/*************************************************************************
 * isoband.cpp
 ************************************************************************/
#include <cmath>
#include <cstring>
#include <map>
#include <vector>

#include "gdal_alg.h"
#include "geos_c.h"
#include "ogr_feature.h"
#include "ogrsf_frmts.h"

/************************************************************************/
/*                         OGRIsobandWriterInfo                         */
/************************************************************************/
typedef struct
{
    OGRLayer *pLayer;
    OGRGeometry *pBoundary;
    GDALProgressFunc pfnProgress;
    void *pProgressArg;
} OGRIsobandWriterInfo;

/************************************************************************/
/*                         GDALIsobandStatistics                        */
/************************************************************************/
typedef struct
{
    int nCount;
    double dfSum;
    double dfMean;
    double dfStdDev;
    double dfMinimum;
    double dfMaximum;
    double dfCoverage;
} GDALIsobandStatistics;

/************************************************************************/
/*                         GDALIsobandClassifierInfo                    */
/************************************************************************/
typedef struct
{
    GDALRasterBand *pBand;
    bool bAllTouched;
    int nInfoToCalculate;
    std::map<int, GDALIsobandStatistics> oStatistics;
    GDALProgressFunc pfnProgress;
    void *pProgressArg;
} GDALIsobandClassifierInfo;

/************************************************************************/
/*                         GDALIsobandWriter                            */
/************************************************************************/
typedef CPLErr (*GDALIsobandWriter)( OGRGeometry *, void * );

/************************************************************************/
/*                         GDALIsobandClassifier                        */
/************************************************************************/
typedef CPLErr (*GDALIsobandClassifier)( int, OGRGeometry **, double *,
                                         void * );

/************************************************************************/
/*                         GDALIsobandGenerator                         */
/************************************************************************/
class GDALIsobandGenerator
{
public:
    GDALIsobandGenerator( GDALIsobandWriter pfnWriter, void *pWriterData,
                          GDALIsobandClassifier pfnClassifier,
                          void *pClassifierData );
    ~GDALIsobandGenerator();
    CPLErr FeedContour( OGRGeometry *pContour );
    CPLErr WriteIsobands();
    CPLErr FeedIsoband( OGRGeometry *pIsoband, int iId );
    CPLErr ClassifyIsobands();

private:
    GDALIsobandWriter pfnWriter;
    void *pWriterData;
    OGRGeometry *pContourUnion;
    
    GDALIsobandClassifier pfnClassifier;
    void *pClassifierData;
    std::map<int, OGRGeometry *> oClonedIsobands;
};

/************************************************************************/
/*                         GDALIsobandGenerator()                       */
/************************************************************************/
GDALIsobandGenerator::GDALIsobandGenerator( GDALIsobandWriter pfnWriterIn,
                                            void *pWriterDataIn,
                                            GDALIsobandClassifier pfnClassifierIn,
                                            void *pClassifierDataIn )
{
    pfnWriter = pfnWriterIn;
    pWriterData = pWriterDataIn;
    pContourUnion = NULL;
    pfnClassifier = pfnClassifierIn;
    pClassifierData = pClassifierDataIn;
}

/************************************************************************/
/*                         ~GDALIsobandGenerator()                      */
/************************************************************************/
GDALIsobandGenerator::~GDALIsobandGenerator()
{
    std::map<int, OGRGeometry *>::iterator it;
    for (it = oClonedIsobands.begin(); it != oClonedIsobands.end(); ++it)
        OGRGeometryFactory::destroyGeometry(it->second);
    OGRGeometryFactory::destroyGeometry(pContourUnion);
}

/************************************************************************/
/*                         FeedContour()                                */
/************************************************************************/
CPLErr GDALIsobandGenerator::FeedContour( OGRGeometry *pContour )
{
    if (pContourUnion)
        pContourUnion = pContourUnion->Union(pContour);
    else
        pContourUnion = pContour->clone();

    if (!pContourUnion)
    {
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not union contour geometry");
        return CE_Failure;
    }
    return CE_None;
}

/************************************************************************/
/*                         WriteIsobands()                              */
/************************************************************************/
CPLErr GDALIsobandGenerator::WriteIsobands()
{
    return pfnWriter(pContourUnion, pWriterData);
}

/************************************************************************/
/*                         FeedIsoband()                                */
/************************************************************************/
CPLErr GDALIsobandGenerator::FeedIsoband( OGRGeometry *pIsoband, int iId )
{
    pIsoband = pIsoband->clone();
    if (!pIsoband)
    {
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not clone isoband geometry");
        return CE_Failure;
    }
    oClonedIsobands[iId] = pIsoband;
    return CE_None;
}

/************************************************************************/
/*                         ClassifyIsobands()                           */
/************************************************************************/
CPLErr GDALIsobandGenerator::ClassifyIsobands()
{
    int nGeomCount = oClonedIsobands.size();
    std::vector<OGRGeometry *> oIsobands;
    std::vector<double> oIds;
    oIsobands.reserve(nGeomCount);
    oIds.reserve(nGeomCount);

    std::map<int, OGRGeometry *>::iterator it;
    for (it = oClonedIsobands.begin(); it != oClonedIsobands.end(); ++it)
    {
        oIsobands.push_back(it->second);
        oIds.push_back(it->first);
    }

    return pfnClassifier(nGeomCount, oIsobands.data(), oIds.data(),
                         pClassifierData);
}

/************************************************************************/
/*                         OGRDefaultIsobandWriter()                    */
/*   Adapted from QGIS: https://qgis.org/api/qgsgeos_8cpp_source.html   */
/************************************************************************/
static CPLErr OGRDefaultIsobandWriter( OGRGeometry *pContours,
                                       void *pInfoIn )
{
    OGRIsobandWriterInfo *pInfo = (OGRIsobandWriterInfo *)pInfoIn;
    void *pScaledProgress;

    if (!pInfo->pfnProgress(0, "", pInfo->pProgressArg))
    {
        CPLError(CE_Failure, CPLE_UserInterrupt, "User terminated");
        return CE_Failure;
    }

    //--------------------------------------------------------------------
    // If no contour was provided, then write the boundary itself
    //--------------------------------------------------------------------
    if (!pContours)
    {
        OGRGeometry *pPolygon = pInfo->pBoundary->clone();
        if (!pPolygon)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "could not clone boundary geometry");
            return CE_Failure;
        }

        OGRFeature *pFeature = OGRFeature::CreateFeature(
            pInfo->pLayer->GetLayerDefn());
        if (!pFeature)
        {
            OGRGeometryFactory::destroyGeometry(pPolygon);
            CPLError(CE_Failure, CPLE_AppDefined,
                     "could not create polygon feature");
            return CE_Failure;
        }

        pFeature->SetGeometryDirectly(pPolygon);

        OGRErr eErr = pInfo->pLayer->CreateFeature(pFeature);
        OGRFeature::DestroyFeature(pFeature);
        if (eErr != OGRERR_NONE)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "could not write polygon to the output layer");
            return CE_Failure;
        }

        if (!pInfo->pfnProgress(1.0, "", pInfo->pProgressArg))
        {
            CPLError(CE_Failure, CPLE_UserInterrupt, "User terminated");
            return CE_Failure;
        }

        return CE_None;
    }

    //--------------------------------------------------------------------
    // Get GEOS representations of the geometries
    //--------------------------------------------------------------------
    GEOSContextHandle_t hGEOSCtxt = OGRGeometry::createGEOSContext();
    GEOSGeom hBoundaryGeosGeom =
        pInfo->pBoundary->exportToGEOS(hGEOSCtxt);
    GEOSGeom hContoursGeosGeom =
        pContours ? pContours->exportToGEOS(hGEOSCtxt) :
                    GEOSGeom_createEmptyLineString_r(hGEOSCtxt);

    if (!hBoundaryGeosGeom || !hContoursGeosGeom)
    {
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not create GEOS geometries");
        return CE_Failure;
    }

    if (!GEOSisValid_r(hGEOSCtxt, hContoursGeosGeom)
        || !GEOSisSimple_r(hGEOSCtxt, hContoursGeosGeom))
    {
        GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
        GEOSGeom_destroy_r(hGEOSCtxt, hContoursGeosGeom);
        OGRGeometry::freeGEOSContext(hGEOSCtxt);
        CPLError(CE_Failure, CPLE_AppDefined,
                 "invalid contour line strings");
        return CE_Failure;
    }

    if (!GEOSIntersects_r(hGEOSCtxt, hContoursGeosGeom,
                          hBoundaryGeosGeom))
    {
        GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
        GEOSGeom_destroy_r(hGEOSCtxt, hContoursGeosGeom);
        OGRGeometry::freeGEOSContext(hGEOSCtxt);
        return CE_None;
    }

    if (!pInfo->pfnProgress(0.1, "", pInfo->pProgressArg))
    {
        CPLError(CE_Failure, CPLE_UserInterrupt, "User terminated");
        return CE_Failure;
    }

    //--------------------------------------------------------------------
    // Union line strings together (node them) and polygonize
    //--------------------------------------------------------------------
    GEOSGeom hBoundaryLineGeosGeom = GEOSBoundary_r(hGEOSCtxt, 
                                                    hBoundaryGeosGeom);
    GEOSGeom hNodedGeometry = GEOSUnion_r(hGEOSCtxt, hContoursGeosGeom,
                                          hBoundaryLineGeosGeom);
    GEOSGeom_destroy_r(hGEOSCtxt, hContoursGeosGeom);
    GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryLineGeosGeom);
    if (!hNodedGeometry)
    {
        GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
        OGRGeometry::freeGEOSContext(hGEOSCtxt);
        CPLError(CE_Failure, CPLE_AppDefined,
                 "could not node geometries");
        return CE_Failure;
    }

    if (!pInfo->pfnProgress(0.2, "", pInfo->pProgressArg))
    {
        CPLError(CE_Failure, CPLE_UserInterrupt, "User terminated");
        return CE_Failure;
    }

    GEOSGeom hPolygons = GEOSPolygonize_r(hGEOSCtxt, &hNodedGeometry, 1);
    GEOSGeom_destroy_r(hGEOSCtxt, hNodedGeometry);
    if (!hPolygons)
    {
        GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
        OGRGeometry::freeGEOSContext(hGEOSCtxt);
        CPLError(CE_Failure, CPLE_AppDefined,
                 "could not polygonize geometries");
        return CE_Failure;
    }

    if (!pInfo->pfnProgress(0.3, "", pInfo->pProgressArg))
    {
        CPLError(CE_Failure, CPLE_UserInterrupt, "User terminated");
        return CE_Failure;
    }

    int nGeometryType = GEOSGeomTypeId_r(hGEOSCtxt, hPolygons);
    if (nGeometryType != GEOS_MULTIPOLYGON
        && nGeometryType != GEOS_GEOMETRYCOLLECTION)
    {
        GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
        OGRGeometry::freeGEOSContext(hGEOSCtxt);
        CPLError(CE_Failure, CPLE_AppDefined,
                 "wrong geometry type after polygonization");
        return CE_Failure;
    }

    //--------------------------------------------------------------------
    // Check which parts fall within the original boundary
    //--------------------------------------------------------------------
    std::vector<GEOSGeom> oGeometriesToAdd;
    int nGeomCount = GEOSGetNumGeometries_r(hGEOSCtxt, hPolygons);
    pScaledProgress = GDALCreateScaledProgress(0.3, 0.7,
                                               pInfo->pfnProgress,
                                               pInfo->pProgressArg);
    for (int iGeom = 0; iGeom < nGeomCount; iGeom++)
    {
        const GEOSGeometry *hPolygon =
            GEOSGetGeometryN_r(hGEOSCtxt, hPolygons, iGeom);
        if (GEOSGeomTypeId_r(hGEOSCtxt, hPolygon) != GEOS_POLYGON)
        {
            GEOSGeom_destroy_r(hGEOSCtxt, hPolygons);
            GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
            OGRGeometry::freeGEOSContext(hGEOSCtxt);
            CPLError(CE_Warning, CPLE_AppDefined,
                     "the generated geometry is not a polygon");
            return CE_Warning;
        }

        GEOSGeom hIntersection = GEOSIntersection_r(hGEOSCtxt,
                                                    hBoundaryGeosGeom,
                                                    hPolygon);
        if (!hIntersection)
        {
            GEOSGeom_destroy_r(hGEOSCtxt, hPolygons);
            GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
            OGRGeometry::freeGEOSContext(hGEOSCtxt);
            CPLError(CE_Warning, CPLE_AppDefined,
                     "could not get intersection");
            return CE_Failure;
        }

        double polygonArea, intersectionArea;
        GEOSArea_r(hGEOSCtxt, hPolygon, &polygonArea);
        GEOSArea_r(hGEOSCtxt, hIntersection, &intersectionArea);
        GEOSGeom_destroy_r(hGEOSCtxt, hIntersection);

        double areaRatio = intersectionArea / polygonArea;
        if (areaRatio > 0.99 && areaRatio < 1.01)
        {
            if (!GEOSisValid_r(hGEOSCtxt, hPolygon))
            {
                GEOSGeom_destroy_r(hGEOSCtxt, hPolygons);
                GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
                OGRGeometry::freeGEOSContext(hGEOSCtxt);
                CPLError(CE_Failure, CPLE_AppDefined,
                         "the generated polygon is invalid");
                return CE_Failure;
            }
            oGeometriesToAdd.push_back(
                GEOSGeom_clone_r(hGEOSCtxt, hPolygon));
        }
        
        if (!GDALScaledProgress(iGeom / (double)nGeomCount, "",
                                pScaledProgress))
        {
            GDALDestroyScaledProgress(pScaledProgress);
            CPLError(CE_Failure, CPLE_UserInterrupt, "User terminated");
            return CE_Failure;
        }
    }
    GEOSGeom_destroy_r(hGEOSCtxt, hPolygons);
    GDALDestroyScaledProgress(pScaledProgress);

    if (oGeometriesToAdd.empty())
    {
        GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
        OGRGeometry::freeGEOSContext(hGEOSCtxt);
        CPLError(CE_Failure, CPLE_AppDefined,
                 "no polygons to add after splitting");
        return CE_Failure;
    }

    //--------------------------------------------------------------------
    // Add the new geometries to the output layer
    //--------------------------------------------------------------------
    pScaledProgress = GDALCreateScaledProgress(0.7, 1.0,
                                               pInfo->pfnProgress,
                                               pInfo->pProgressArg);
    for (std::size_t iGeom = 0; iGeom < oGeometriesToAdd.size(); iGeom++)
    {
        OGRGeometry *pPolygon =
            OGRGeometryFactory::createFromGEOS(hGEOSCtxt,
                                               oGeometriesToAdd[iGeom]);
        GEOSGeom_destroy_r(hGEOSCtxt, oGeometriesToAdd[iGeom]);

        if (!pPolygon)
        {
            GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
            OGRGeometry::freeGEOSContext(hGEOSCtxt);
            CPLError(CE_Failure, CPLE_AppDefined,
                     "could not create OGR geometry from GEOS geometry");
            return CE_Failure;
        }
        pPolygon->assignSpatialReference(
            pInfo->pBoundary->getSpatialReference());

        OGRFeature *pFeature = OGRFeature::CreateFeature(
            pInfo->pLayer->GetLayerDefn());
        if (!pFeature)
        {
            OGRGeometryFactory::destroyGeometry(pPolygon);
            GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
            OGRGeometry::freeGEOSContext(hGEOSCtxt);
            CPLError(CE_Failure, CPLE_AppDefined,
                     "could not create polygon feature");
            return CE_Failure;
        }

        pFeature->SetGeometryDirectly(pPolygon);

        OGRErr eErr = pInfo->pLayer->CreateFeature(pFeature);
        OGRFeature::DestroyFeature(pFeature);
        if (eErr != OGRERR_NONE)
        {
            GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
            OGRGeometry::freeGEOSContext(hGEOSCtxt);
            CPLError(CE_Failure, CPLE_AppDefined,
                     "could not write polygon to the output layer");
            return CE_Failure;
        }
        
        if (!GDALScaledProgress(iGeom / (double)nGeomCount, "",
                                pScaledProgress))
        {
            GDALDestroyScaledProgress(pScaledProgress);
            CPLError(CE_Failure, CPLE_UserInterrupt, "User terminated");
            return CE_Failure;
        }
    }
    GDALDestroyScaledProgress(pScaledProgress);

    GEOSGeom_destroy_r(hGEOSCtxt, hBoundaryGeosGeom);
    OGRGeometry::freeGEOSContext(hGEOSCtxt);

    return CE_None;
}

/************************************************************************/
/*                         AccumulateStatistics()                       */
/************************************************************************/
static void AccumulateStatistics( int iIsobandId, double dfPixelValue,
                                  GDALIsobandClassifierInfo *pInfo )
{
    if (iIsobandId == OGRNullFID)
        return;

    //--------------------------------------------------------------------
    // Initialize statistics
    //--------------------------------------------------------------------
    std::map<int, GDALIsobandStatistics>::iterator it = 
        pInfo->oStatistics.find(iIsobandId);
    if (it == pInfo->oStatistics.end())
    {
        GDALIsobandStatistics oStats =
            {0, 0, NAN, NAN, INFINITY, -INFINITY, NAN};
        it = pInfo->oStatistics.insert(
            std::make_pair(iIsobandId, oStats)).first;
    }

    //--------------------------------------------------------------------
    // Check calculate standard deviation
    //--------------------------------------------------------------------
    bool bCalculateStdDev =
        (pInfo->nInfoToCalculate & GIST_StdDev) == GIST_StdDev;
    if (bCalculateStdDev && !std::isnan(it->second.dfMean))
    {
        double dfDifference = dfPixelValue - it->second.dfMean;
        it->second.dfSum += dfDifference * dfDifference;
        return;
    }

    //--------------------------------------------------------------------
    // Check calculate mean
    //--------------------------------------------------------------------
    bool bCalculateMean =
        (pInfo->nInfoToCalculate & GIST_Mean) == GIST_Mean;
    if (bCalculateMean || bCalculateStdDev)
    {
        ++it->second.nCount;
        it->second.dfSum += dfPixelValue;
    }

    //--------------------------------------------------------------------
    // Check calculate minimum and maximum
    //--------------------------------------------------------------------
    bool bCalculateMinimum =
        (pInfo->nInfoToCalculate & GIST_Minimum) == GIST_Minimum;
    if (bCalculateMinimum && dfPixelValue < it->second.dfMinimum)
        it->second.dfMinimum = dfPixelValue;

    bool bCalculateMaximum =
        (pInfo->nInfoToCalculate & GIST_Maximum) == GIST_Maximum;
    if (bCalculateMaximum && dfPixelValue > it->second.dfMaximum)
        it->second.dfMaximum = dfPixelValue;
}

/************************************************************************/
/*                         IterateOverPixelValues()                     */
/************************************************************************/
template <class T>
CPLErr IterateOverPixelValues( int *pabyRasterized, T *pabyChunkBuffer,
                               int nXSize, int nYSize, int nXBlockSize,
                               int nYBlockSize,
                               GDALIsobandClassifierInfo *pInfo )
{
    CPLErr eErr = CE_None;
    int nXBlockCount = (nXSize + nXBlockSize - 1) / nXBlockSize;
    int nYBlockCount = (nYSize + nYBlockSize - 1) / nYBlockSize;

    //--------------------------------------------------------------------
    // Loop over raster blocks
    //--------------------------------------------------------------------
    for (int i = 0; i < nXBlockCount && eErr == CE_None; ++i)
    {
        for (int j = 0; j < nYBlockCount && eErr == CE_None; ++j)
        {
            eErr = pInfo->pBand->ReadBlock(i, j, pabyChunkBuffer);
            if (eErr != CE_None)
                break;

            //------------------------------------------------------------
            // Check valid pixels within this block
            //------------------------------------------------------------
            int nXValid = (i + 1) * nXBlockSize > nXSize
                ? nXSize - i * nXBlockSize : nXBlockSize;
            int nYValid = (j + 1) * nYBlockSize > nYSize
                ? nYSize - j * nYBlockSize : nYBlockSize;
            for (int k = 0; k < nXValid; ++k)
            {
                for (int l = 0; l < nYValid; ++l)
                {
                    int m = i * nXBlockSize +
                        (j * nYBlockSize + l) * nXSize + k;
                    int n = l * nXBlockSize + k;
                    AccumulateStatistics(pabyRasterized[m],
                                         pabyChunkBuffer[n], pInfo);
                }
            }
        }
    }
    return eErr;
}

/************************************************************************/
/*              IterateOverPixelValuesDependingOnDataType()             */
/************************************************************************/
static CPLErr IterateOverPixelValuesDependingOnDataType(
    int *pabyRasterized, GByte *pabyChunkBuffer, int nXSize, int nYSize,
    int nXBlockSize, int nYBlockSize, GDALIsobandClassifierInfo *pInfo )
{
    switch (pInfo->pBand->GetRasterDataType())
    {
    case GDT_Byte:
        return IterateOverPixelValues(
            pabyRasterized, pabyChunkBuffer, nXSize, nYSize, nXBlockSize,
            nYBlockSize, pInfo);
    case GDT_UInt16:
        return IterateOverPixelValues(
            pabyRasterized, (GUInt16 *)pabyChunkBuffer, nXSize, nYSize,
            nXBlockSize, nYBlockSize, pInfo);
    case GDT_Int16:
        return IterateOverPixelValues(
            pabyRasterized, (GInt16 *)pabyChunkBuffer, nXSize, nYSize,
            nXBlockSize, nYBlockSize, pInfo);
    case GDT_UInt32:
        return IterateOverPixelValues(
            pabyRasterized, (GUInt32 *)pabyChunkBuffer, nXSize, nYSize,
            nXBlockSize, nYBlockSize, pInfo);
    case GDT_Int32:
        return IterateOverPixelValues(
            pabyRasterized, (GInt32 *)pabyChunkBuffer, nXSize, nYSize,
            nXBlockSize, nYBlockSize, pInfo);
    case GDT_Float32:
        return IterateOverPixelValues(
            pabyRasterized, (float *)pabyChunkBuffer, nXSize, nYSize,
            nXBlockSize, nYBlockSize, pInfo);
    case GDT_Float64:
        return IterateOverPixelValues(
            pabyRasterized, (double *)pabyChunkBuffer, nXSize, nYSize,
            nXBlockSize, nYBlockSize, pInfo);
    default:
        CPLError(CE_Failure, CPLE_AppDefined,
                 "raster data type not allowed");
        return CE_Failure;
    }
}

/************************************************************************/
/*              AllocateChunkBufferDependingOnDataType()                */
/************************************************************************/
static void *AllocateChunkBufferDependingOnDataType( int nXBlockSize,
                                                     int nYBlockSize,
                                                     GDALIsobandClassifierInfo *pInfo )
{
    switch (pInfo->pBand->GetRasterDataType())
    {
    case GDT_Byte:
        return CPLCalloc(nXBlockSize * nYBlockSize, 1);
    case GDT_UInt16:
        return CPLCalloc(nXBlockSize * nYBlockSize, sizeof(GUInt16));
    case GDT_Int16:
        return CPLCalloc(nXBlockSize * nYBlockSize, sizeof(GInt16));
    case GDT_UInt32:
        return CPLCalloc(nXBlockSize * nYBlockSize, sizeof(GUInt32));
    case GDT_Int32:
        return CPLCalloc(nXBlockSize * nYBlockSize, sizeof(GInt32));
    case GDT_Float32:
        return CPLCalloc(nXBlockSize * nYBlockSize, sizeof(float));
    case GDT_Float64:
        return CPLCalloc(nXBlockSize * nYBlockSize, sizeof(double));
    default:
        return NULL;
    }
}

/************************************************************************/
/*                     OGRDefaultIsobandClassifier()                    */
/************************************************************************/
static CPLErr OGRDefaultIsobandClassifier( int nGeomCount,
                                           OGRGeometry **papoIsobands,
                                           double *padfClasses,
                                           void *pInfoIn )
{
    GDALIsobandClassifierInfo *pInfo =
        (GDALIsobandClassifierInfo *)pInfoIn;

    //--------------------------------------------------------------------
    // Use memory driver to create a dataset to hold the rasterized
    // geometries
    //--------------------------------------------------------------------
    GDALDriverH hMemoryDriver = GDALGetDriverByName("MEM");
    if (!hMemoryDriver)
    {
        CPLError(CE_Failure, CPLE_ObjectNull, "could not get MEM driver");
        return CE_Failure;
    }

    //--------------------------------------------------------------------
    // Allocate buffer to use as the underlying buffer and initialize it
    //--------------------------------------------------------------------
    int nXSize = pInfo->pBand->GetXSize();
    int nYSize = pInfo->pBand->GetYSize();
    int *pabyRasterized = (int *)CPLCalloc(nXSize * nYSize, sizeof(int));
    if (!pabyRasterized)
    {
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not create buffer for rasterized geometries");
        return CE_Failure;
    }

    for (int i = 0; i < nXSize * nYSize; ++i)
        pabyRasterized[i] = OGRNullFID;

    char szDataPointer[100];
    char *apszOptions[] = {szDataPointer, NULL};
    std::memset(szDataPointer, 0, sizeof(szDataPointer));
    snprintf(szDataPointer, sizeof(szDataPointer), "DATAPOINTER=");
    CPLPrintPointer(szDataPointer + strlen(szDataPointer), pabyRasterized,
                    static_cast<int>(sizeof(szDataPointer) -
                        std::strlen(szDataPointer)));

    //--------------------------------------------------------------------
    // Create the memory dataset as 32-bit integer with one band
    //--------------------------------------------------------------------
    GDALDataset *pMemoryDataset = (GDALDataset *)
        GDALCreate(hMemoryDriver, "", nXSize, nYSize, 0, GDT_Int32, NULL);
    if (!pMemoryDataset)
    {
        CPLFree(pabyRasterized);
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not create memory dataset");
        return CE_Failure;
    }

    CPLErr eErr = pMemoryDataset->AddBand(GDT_Int32, apszOptions);
    if (eErr != CE_None)
    {
        GDALClose(pMemoryDataset);
        CPLFree(pabyRasterized);
        return eErr;
    }

    //--------------------------------------------------------------------
    // Assign geo transform
    //--------------------------------------------------------------------
    GDALDataset *pDataset = pInfo->pBand->GetDataset();
    if (!pDataset)
    {
        GDALClose(pMemoryDataset);
        CPLFree(pabyRasterized);
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not get band dataset");
        return CE_Failure;
    }

    double adfGeoTransform[6];
    pDataset->GetGeoTransform(adfGeoTransform);
    pMemoryDataset->SetGeoTransform(adfGeoTransform);

    //--------------------------------------------------------------------
    // Rasterize isobands with burn values corresponding to their indices
    //--------------------------------------------------------------------
    int nTargetBand = 1;
    char **papszRasterizeOptions = NULL;
    if (pInfo->bAllTouched)
        papszRasterizeOptions = CSLSetNameValue(papszRasterizeOptions,
                                                "ALL_TOUCHED", "TRUE");

    OGRGeometryH *pahIsobands = (OGRGeometryH *)papoIsobands;
    eErr = GDALRasterizeGeometries(pMemoryDataset, 1, &nTargetBand,
                                   nGeomCount, pahIsobands, NULL, NULL,
                                   padfClasses, papszRasterizeOptions,
                                   pInfo->pfnProgress,
                                   pInfo->pProgressArg);
    CSLDestroy(papszRasterizeOptions);
    GDALClose(pMemoryDataset);

    if (eErr != CE_None)
    {
        CPLFree(pabyRasterized);
        return eErr;
    }

    //--------------------------------------------------------------------
    // Allocate buffer for reading the input raster
    //--------------------------------------------------------------------
    int nXBlockSize, nYBlockSize;
    pInfo->pBand->GetBlockSize(&nXBlockSize, &nYBlockSize);
    GByte *pabyChunkBuffer = (GByte *)
        AllocateChunkBufferDependingOnDataType(nXBlockSize, nYBlockSize,
                                               pInfo);
    if (!pabyChunkBuffer)
    {
        CPLFree(pabyRasterized);
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not create buffer for reading input raster");
        return CE_Failure;
    }

    //--------------------------------------------------------------------
    // Compute the required statistics (first pass)
    //--------------------------------------------------------------------
    eErr = IterateOverPixelValuesDependingOnDataType(
        pabyRasterized, pabyChunkBuffer, nXSize, nYSize, nXBlockSize,
        nYBlockSize, pInfo);
    if (eErr != CE_None)
    {
        CPLFree(pabyChunkBuffer);
        CPLFree(pabyRasterized);
        return eErr;
    }

    std::map<int, GDALIsobandStatistics>::iterator it;
    for (it = pInfo->oStatistics.begin();
        it != pInfo->oStatistics.end(); ++it)
    {
        it->second.dfMean = it->second.dfSum / (double)it->second.nCount;
        it->second.dfSum = 0;
    }

    //--------------------------------------------------------------------
    // Compute the required statistics (second pass)
    //--------------------------------------------------------------------
    if ((pInfo->nInfoToCalculate & GIST_StdDev) == GIST_StdDev)
    {
        eErr = IterateOverPixelValuesDependingOnDataType(
            pabyRasterized, pabyChunkBuffer, nXSize, nYSize, nXBlockSize,
            nYBlockSize, pInfo);

        if (eErr == CE_None)
        {
            for (it = pInfo->oStatistics.begin();
                it != pInfo->oStatistics.end(); ++it)
                it->second.dfStdDev =
                    std::sqrt(it->second.dfSum) /
                        (double)it->second.nCount;
        }
    }
    CPLFree(pabyChunkBuffer);
    CPLFree(pabyRasterized);

    return eErr;
}

/************************************************************************/
/*                          GetBoundaryFromLayer()                      */
/************************************************************************/
static CPLErr GetBoundaryFromLayer( OGRLayer *pLayer,
                                    OGRPolygon **papoResult )
{
    pLayer->ResetReading();

    OGRFeature *pBoundaryFeature = pLayer->GetNextFeature();
    if (!pBoundaryFeature)
    {
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "no features found in the boundary layer");
        return CE_Failure;
    }

    OGRGeometry *pBoundaryGeometry = pBoundaryFeature->GetGeometryRef();
    if (!pBoundaryGeometry)
    {
        OGRFeature::DestroyFeature(pBoundaryFeature);
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "the boundary feature is empty");
        return CE_Failure;
    }

    if (wkbFlatten(pBoundaryGeometry->getGeometryType()) != wkbPolygon)
    {
        OGRFeature::DestroyFeature(pBoundaryFeature);
        CPLError(CE_Failure, CPLE_AppDefined,
                 "the boundary geometry must be a polygon");
        return CE_Failure;
    }

    if (!pBoundaryGeometry->IsValid())
    {
        OGRFeature::DestroyFeature(pBoundaryFeature);
        CPLError(CE_Failure, CPLE_AppDefined,
                 "invalid boundary polygon");
        return CE_Failure;
    }

    pBoundaryGeometry = pBoundaryGeometry->clone();
    OGRFeature::DestroyFeature(pBoundaryFeature);
    if (!papoResult)
    {
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not clone the boundary polygon");
        return CE_Failure;
    }

    *papoResult = (OGRPolygon *)pBoundaryGeometry;
    return CE_None;
}

/************************************************************************/
/*                    CreateBoundaryFromRasterExtent()                  */
/************************************************************************/
static CPLErr CreateBoundaryFromRasterExtent( GDALRasterBand *pBand,
                                              OGRPolygon **papoResult )
{
    *papoResult = NULL;

    //--------------------------------------------------------------------
    // Get raster geotransform parameters
    //--------------------------------------------------------------------
    GDALDataset *pDataset = pBand->GetDataset();
    double adfGeoTransform[6];
    pDataset->GetGeoTransform(adfGeoTransform);

    double minX = adfGeoTransform[0];
    double maxX = pDataset->GetRasterXSize() * adfGeoTransform[1] + minX;
    double maxY = adfGeoTransform[3];
    double minY = pDataset->GetRasterYSize() * adfGeoTransform[5] + maxY;

    //--------------------------------------------------------------------
    // Create and fill a ring to add to the polygon
    //--------------------------------------------------------------------
    OGRLinearRing *pRing = (OGRLinearRing *)
        OGRGeometryFactory::createGeometry(wkbLinearRing);
    if (!pRing)
    {
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not create line string geometry");
        return CE_Failure;
    }

    pRing->addPoint(minX, maxY);
    pRing->addPoint(maxX, maxY);
    pRing->addPoint(maxX, minY);
    pRing->addPoint(minX, minY);
    pRing->addPoint(minX, maxY);

    //--------------------------------------------------------------------
    // Create and fill the resulting polygon
    //--------------------------------------------------------------------
    OGRPolygon *pPolygon = (OGRPolygon *)
        OGRGeometryFactory::createGeometry(wkbPolygon);
    if (!pPolygon)
    {
        CPLError(CE_Failure, CPLE_ObjectNull,
                 "could not create polygon geometry");
        return CE_Failure;
    }

    OGRErr eErr = pPolygon->addRing(pRing);
    OGRGeometryFactory::destroyGeometry(pRing);
    if (eErr != OGRERR_NONE)
    {
        OGRGeometryFactory::destroyGeometry(pPolygon);
        CPLError(CE_Failure, CPLE_AppDefined,
                 "could not add ring to polygon");
        return CE_Failure;
    }

    *papoResult = pPolygon;
    return CE_None;
}

/************************************************************************/
/*                             CreateField()                            */
/************************************************************************/
static CPLErr CreateField( OGRLayer *pLayer, const char *pszName,
                           OGRFieldType nType )
{
    OGRFieldDefn oDefinition(pszName, nType);
    OGRErr eErr = pLayer->CreateField(&oDefinition);
    if (eErr != OGRERR_NONE)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "could not create field in output layer");
        return CE_Failure;
    }
    return CE_None;
}

/************************************************************************/
/*                            CreateFields()                            */
/************************************************************************/
static CPLErr CreateFields( int nInfoToCalculate, OGRLayer *pLayer,
                            char **papszOptions )
{
    CPLErr eErr = CE_None;
    if ((nInfoToCalculate & GIST_SampleCount) == GIST_SampleCount)
        eErr = CreateField(
            pLayer, CSLFetchNameValueDef(
                papszOptions, "SAMPLE_COUNT_FIELD_NAME", "count"),
                OFTInteger);
    if (eErr == CE_None && (nInfoToCalculate & GIST_Mean) == GIST_Mean)
        eErr = CreateField(
            pLayer, CSLFetchNameValueDef(
                papszOptions, "MEAN_FIELD_NAME", "mean"), OFTReal);
    if (eErr == CE_None
        && (nInfoToCalculate & GIST_StdDev) == GIST_StdDev)
        eErr = CreateField(
            pLayer, CSLFetchNameValueDef(
                papszOptions, "STDDEV_FIELD_NAME", "stddev"), OFTReal);
    if (eErr == CE_None
        && (nInfoToCalculate & GIST_Minimum) == GIST_Minimum)
        eErr = CreateField(
            pLayer, CSLFetchNameValueDef(
                papszOptions, "MINIMUM_FIELD_NAME", "minimum"), OFTReal);
    if (eErr == CE_None
        && (nInfoToCalculate & GIST_Maximum) == GIST_Maximum)
        eErr = CreateField(
            pLayer, CSLFetchNameValueDef(
                papszOptions, "MAXIMUM_FIELD_NAME", "maximum"), OFTReal);
    if (eErr == CE_None
        && (nInfoToCalculate & GIST_Coverage) == GIST_Coverage)
        eErr = CreateField(
            pLayer, CSLFetchNameValueDef(
                papszOptions, "COVERAGE_FIELD_NAME", "coverage"),
                OFTReal);
    return eErr;
}

/************************************************************************/
/*                            WriteFields()                             */
/************************************************************************/
static void WriteFields( int nInfoToCalculate, OGRFeature *pFeature,
                         GDALIsobandStatistics *pStats,
                         char **papszOptions )
{
    if ((nInfoToCalculate & GIST_SampleCount) == GIST_SampleCount)
        pFeature->SetField(
            CSLFetchNameValueDef(
                papszOptions, "SAMPLE_COUNT_FIELD_NAME", "count"),
                pStats->nCount);
    if ((nInfoToCalculate & GIST_Mean) == GIST_Mean)
        pFeature->SetField(
            CSLFetchNameValueDef(
                papszOptions, "MEAN_FIELD_NAME", "mean"), pStats->dfMean);
    if ((nInfoToCalculate & GIST_StdDev) == GIST_StdDev)
        pFeature->SetField(
            CSLFetchNameValueDef(
                papszOptions, "STDDEV_FIELD_NAME", "stddev"),
                pStats->dfStdDev);
    if ((nInfoToCalculate & GIST_Minimum) == GIST_Minimum)
        pFeature->SetField(
            CSLFetchNameValueDef(
                papszOptions, "MINIMUM_FIELD_NAME", "minimum"),
                pStats->dfMinimum);
    if ((nInfoToCalculate & GIST_Maximum) == GIST_Maximum)
        pFeature->SetField(
            CSLFetchNameValueDef(
                papszOptions, "MAXIMUM_FIELD_NAME", "maximum"),
                pStats->dfMaximum);
    if ((nInfoToCalculate & GIST_Coverage) == GIST_Coverage)
        pFeature->SetField(
            CSLFetchNameValueDef(
                papszOptions, "COVERAGE_FIELD_NAME", "coverage"),
                pStats->dfCoverage);
}

/************************************************************************/
/*                         GDALIsobandGenerate()                        */
/************************************************************************/
CPLErr GDALIsobandGenerate( GDALRasterBandH hRasterBand,
                            void *hContourLayer, void *hBoundaryLayer,
                            void *hIsobandLayer, int nInfoToCalculate,
                            int bAllTouched, char **papszOptions,
                            GDALProgressFunc pfnProgress,
                            void *pProgressArg)
{
    CPLErr eErr;

    //--------------------------------------------------------------------
    // Get polygon geometry from boundary layer or create one based on
    // the raster extent
    //--------------------------------------------------------------------
    OGRPolygon *pBoundaryPolygon;
    if (hBoundaryLayer)
        eErr = GetBoundaryFromLayer((OGRLayer *)hBoundaryLayer,
                                    &pBoundaryPolygon);
    else
        eErr = CreateBoundaryFromRasterExtent(
            (GDALRasterBand *)hRasterBand, &pBoundaryPolygon);

    if (eErr != CE_None)
        return eErr;

    //--------------------------------------------------------------------
    // Create required fields
    //--------------------------------------------------------------------
    eErr = CreateFields(nInfoToCalculate, (OGRLayer *)hIsobandLayer,
                        papszOptions);
    if (eErr != CE_None)
    {
        OGRGeometryFactory::destroyGeometry(pBoundaryPolygon);
        return eErr;
    }

    //--------------------------------------------------------------------
    // Setup isoband generator
    //--------------------------------------------------------------------
    OGRIsobandWriterInfo oIWI;
    oIWI.pLayer = (OGRLayer *)hIsobandLayer;
    oIWI.pBoundary = pBoundaryPolygon;

    GDALIsobandClassifierInfo oICI;
    oICI.pBand = (GDALRasterBand *)hRasterBand;
    oICI.bAllTouched = bAllTouched;
    oICI.nInfoToCalculate = nInfoToCalculate;

    GDALIsobandGenerator oIG(OGRDefaultIsobandWriter, &oIWI,
                             OGRDefaultIsobandClassifier, &oICI);

    //--------------------------------------------------------------------
    // Feed the contours into the generator
    //--------------------------------------------------------------------
    OGRLayer *pContourLayer = (OGRLayer *)hContourLayer;
    std::size_t iContour = 0;
    std::size_t nContours = pContourLayer->GetFeatureCount();

    pContourLayer->ResetReading();
    OGRFeature *pContourFeature;

    void *pScaledProgress = GDALCreateScaledProgress(0, 0.2, pfnProgress,
                                                     pProgressArg);
    while ((pContourFeature = pContourLayer->GetNextFeature())
           && eErr == CE_None)
    {
        eErr = oIG.FeedContour(pContourFeature->GetGeometryRef());
        OGRFeature::DestroyFeature(pContourFeature);

        if (eErr == CE_None
            && !GDALScaledProgress(++iContour / (double)nContours, "",
                                   pScaledProgress))
        {
            CPLError(CE_Failure, CPLE_UserInterrupt, "User terminated");
            eErr = CE_Failure;
        }
    }
    GDALDestroyScaledProgress(pScaledProgress);

    //--------------------------------------------------------------------
    // Split the boundary and write the isobands
    //--------------------------------------------------------------------
    double dfProgressEnd = (nInfoToCalculate != GIST_None) ? 0.6 : 1.0;
    oIWI.pfnProgress = GDALScaledProgress;
    oIWI.pProgressArg = GDALCreateScaledProgress(0.2, dfProgressEnd,
                                                 pfnProgress,
                                                 pProgressArg);
    eErr = oIG.WriteIsobands();
    GDALDestroyScaledProgress(oIWI.pProgressArg);

    if (eErr == CE_None && nInfoToCalculate != GIST_None)
    {
        //----------------------------------------------------------------
        // Feed the isobands into the generator (for classification)
        //----------------------------------------------------------------
        OGRLayer *pIsobandLayer = (OGRLayer *)hIsobandLayer;
        std::size_t iIsoband = 0;
        std::size_t nIsobands = pIsobandLayer->GetFeatureCount();

        pIsobandLayer->ResetReading();
        OGRFeature *pIsobandFeature;

        pScaledProgress = GDALCreateScaledProgress(0.6, 0.7, pfnProgress,
                                                   pProgressArg);
        while ((pIsobandFeature = pIsobandLayer->GetNextFeature())
            && eErr == CE_None)
        {
            eErr = oIG.FeedIsoband(pIsobandFeature->GetGeometryRef(),
                                   pIsobandFeature->GetFID());
            OGRFeature::DestroyFeature(pIsobandFeature);

            if (eErr == CE_None
                && !GDALScaledProgress(++iIsoband / (double)nIsobands, "",
                                       pScaledProgress))
            {
                CPLError(CE_Failure, CPLE_UserInterrupt,
                         "User terminated");
                eErr = CE_Failure;
            }
        }
        GDALDestroyScaledProgress(pScaledProgress);

        if (eErr != CE_None)
        {
            OGRGeometryFactory::destroyGeometry(pBoundaryPolygon);
            return eErr;
        }

        //----------------------------------------------------------------
        // Classify isobands and compute statistics
        //----------------------------------------------------------------
        oICI.pfnProgress = GDALScaledProgress;
        oICI.pProgressArg = GDALCreateScaledProgress(0.7, 0.9,
                                                     pfnProgress,
                                                     pProgressArg);
        oIG.ClassifyIsobands();
        GDALDestroyScaledProgress(oICI.pProgressArg);

        double dfTotalArea = 0;
        if ((nInfoToCalculate & GIST_Coverage) == GIST_Coverage)
            dfTotalArea = pBoundaryPolygon->get_Area();

        //----------------------------------------------------------------
        // Write the statistics into appropriate fields
        //----------------------------------------------------------------
        iIsoband = 0;
        nIsobands = oICI.oStatistics.size();
        std::map<int, GDALIsobandStatistics>::iterator it;

        pScaledProgress = GDALCreateScaledProgress(0.9, 1.0, pfnProgress,
                                                   pProgressArg);
        for (it = oICI.oStatistics.begin();
            it != oICI.oStatistics.end(); ++it)
        {
            OGRFeature *poFeature = pIsobandLayer->GetFeature(it->first);
            if (!poFeature)
            {
                CPLError(CE_Failure, CPLE_ObjectNull,
                         "could not find the required feature");
                eErr = CE_Failure;
                break;
            }

            OGRPolygon *pIsoband =
                (OGRPolygon *)poFeature->GetGeometryRef();
            if (!pIsoband)
            {
                OGRFeature::DestroyFeature(poFeature);
                CPLError(CE_Failure, CPLE_ObjectNull,
                         "could not get the geometry of this feature");
                eErr = CE_Failure;
                break;
            }

            if ((nInfoToCalculate & GIST_Coverage) == GIST_Coverage)
                it->second.dfCoverage =
                    pIsoband->get_Area() / dfTotalArea;

            WriteFields(nInfoToCalculate, poFeature, &it->second,
                        papszOptions);

            OGRErr eOgrErr = pIsobandLayer->SetFeature(poFeature);
            OGRFeature::DestroyFeature(poFeature);
            if (eOgrErr != OGRERR_NONE)
            {
                CPLError(CE_Failure, CPLE_AppDefined,
                         "could not write feature to output layer");
                eErr = CE_Failure;
                break;
            }

            if (!GDALScaledProgress(++iIsoband / (double)nIsobands, "",
                                    pScaledProgress))
            {
                CPLError(CE_Failure, CPLE_UserInterrupt,
                         "User terminated");
                eErr = CE_Failure;
                break;
            }
        }
        GDALDestroyScaledProgress(pScaledProgress);
    }
    OGRGeometryFactory::destroyGeometry(pBoundaryPolygon);

    return eErr;
}
