// summer_pra_4_1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <string>
#include "gdal201/gdal/include/ogrsf_frmts.h"
#include "gdal201/gdal/include/gdal_priv.h"
#pragma comment(lib, "gdal_i.lib")
using namespace std;

typedef struct OutRing {
    vector<double> Ringx;
    vector<double> Ringy;
}OutRing;
typedef struct Point {
    double pointx;
    double pointy;
}Point;

bool isequal(Point a, Point b);


int main() {
    //cout << "hello";

    GDALAllRegister();
    CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
    CPLSetConfigOption("SHAPE_ENCODING", "");

    GDALDataset* pDataset = NULL;
    pDataset = (GDALDataset*)GDALOpenEx("data/Demo_Polygon.shp", GDAL_OF_VECTOR, NULL, NULL, NULL);
    //创建句柄
    int iLayerNum = pDataset->GetLayerCount();
    //获取层数
    OGRLayer* pLayer = pDataset->GetLayer(0);
    //读取第一层
    OGREnvelope Envelope;
    pDataset->GetLayer(0)->GetExtent(&Envelope);
    //获取包围盒
    string strLayerName = pLayer->GetName();
    //获取图层名
    OGRSpatialReference* SpatialReference = pLayer->GetSpatialRef();
    //获取空间参考
    OGRwkbGeometryType GeometryType = pLayer->GetGeomType();
    //获取图元类型
    int64_t iFeatureNum = pLayer->GetFeatureCount();
    //获取要素数量
    cout << iFeatureNum << endl;//调试
    OGRFeatureDefn* pAttribute = pLayer->GetLayerDefn();
    //
    int iFieldNum = pAttribute->GetFieldCount();
    //
    OGRFieldType iFieldType;
    string iFieldName;
    vector<OGRFieldType> types;
    vector<string> names;


    OGRFieldDefn* pField = NULL;
    for (int i = 0; i < iFieldNum; i++) {
        pField = pAttribute->GetFieldDefn(i);
        iFieldType = pField->GetType();
        iFieldName = pField->GetNameRef();
        types.push_back(iFieldType);
        names.push_back(iFieldName);
        cout << names[i] << "\t" << types[i] << endl;
    }

    vector<vector<string>> FieldsValues_list;
    OGRFeature* pFeature;
    vector<string> FieldsValue;
    vector<OutRing> outrings;
    vector<vector<OutRing>> inholes;

    for (int j = 0; j < iFeatureNum; j++) {

        //field readin
        pFeature = pLayer->GetFeature(j);

        FieldsValue.clear();

        for (int i = 0; i < iFieldNum; i++) {
            string strValue = pFeature->GetFieldAsString(i);
            FieldsValue.push_back(strValue);
            cout << FieldsValue[i] << endl;
        }
        FieldsValues_list.push_back(FieldsValue);


        //geometry readin
        OGRGeometry* pGeometry = pFeature->GetGeometryRef();
        //OGRwkbGeometryType type = wkbFlatten(pGeometry->getGeometryType());//判断是什么图元类型，本次实习全是多边形所以直接不用这句话
        OGRPolygon* pOGRPolygon = (OGRPolygon*)pGeometry;

        OGRLinearRing* pRing = pOGRPolygon->getExteriorRing();
        int Ringpointnum = pRing->getNumPoints();
        OutRing rings;
        vector<double> ringx;
        vector<double> ringy;
        for (int i = 0; i < Ringpointnum; i++) {
            ringx.push_back(pRing->getX(i));
            ringy.push_back(pRing->getY(i));
        }
        rings.Ringx = ringx; rings.Ringy = ringy;
        outrings.push_back(rings);
        //readin done (rings)
        for (int i = 0; i < Ringpointnum; i++) {
            cout << rings.Ringx[i] << "\t" << rings.Ringy[i] << endl;
        }//hint output

        cout << endl;
        vector<double> holepx;
        vector<double> holepy;
        int Holenum = pOGRPolygon->getNumInteriorRings();
        int Holepointnum;
        vector<OutRing> holes;
        OutRing holetmp;
        OGRLinearRing* pHole;
        for (int i = 0; i < Holenum; i++) {
            pHole = pOGRPolygon->getInteriorRing(i);
            Holepointnum = pHole->getNumPoints();
            for (int k = 0; k < Holepointnum; k++) {
                holepx.push_back(pHole->getX(k));
                holepy.push_back(pHole->getY(k));
            }
            holetmp.Ringx = holepx;
            holetmp.Ringy = holepy;
            holes.push_back(holetmp);
            holepx.clear();
            holepy.clear();
        }
        inholes.push_back(holes);
        //readin done (holes)
        for (int i = 0; i < Holenum; i++) {
            pHole = pOGRPolygon->getInteriorRing(i);
            Holepointnum = pHole->getNumPoints();
            for (int k = 0; k < Holepointnum; k++) {
                cout << holes[i].Ringx[k] << "\t" << holes[i].Ringy[k] << endl;
            }
            cout << i << endl;
        }//hint output
    }

    int region1, region2;
    cin >> region1 >> region2;
    Point pointtmp1, pointtmp2;
    vector<int> repeatp1, repeatp2;
    int resize = 0;
    bool flag = false;
    vector<int> respolygon;
    int regionsize = 0;
    //becauseof the 0 and the size of outrings, the mark is begin from 1
    //search in outrings we always begin with 1
    //judge if the two region is near and record the same points
    for (int i = 1; i < outrings[region1].Ringx.size(); i++) {
        pointtmp1.pointx = outrings[region1].Ringx[i];
        pointtmp1.pointy = outrings[region1].Ringy[i];
        for (int j = 1; j < outrings[region2].Ringy.size(); j++) {
            pointtmp2.pointx = outrings[region2].Ringx[j];
            pointtmp2.pointy = outrings[region2].Ringy[j];
            if (isequal(pointtmp1, pointtmp2)) {
                flag = true;
                resize++;
                repeatp1.push_back(i);
                repeatp2.push_back(j);
                break;
            }
        }
    }
    //start to search two ring, negative number mark means region2
    int beginnum = 0;
    int nextnum = 0;
    if (flag) {
        //find a not same point in region1
        bool mark = true;
        for (int i = 1; i < outrings[region1].Ringx.size(); i++) {
            beginnum = 0;
            for (int j = 0; j < resize; j++) {
                if (i == repeatp1[j]) {
                    beginnum = i;
                    break;
                }
            }
            if (beginnum == 0) {
                beginnum = i;
                respolygon.push_back(i);
                regionsize++;
                nextnum = beginnum;
                break;
            }
        }
        //start push pointlist in resregion
        //this algorithm is unfittable for all of the situations like "inneroutring", "all samepoint region" etc.
        while (true) {
            if (mark) {
                nextnum++;
                if (nextnum == outrings[region1].Ringx.size()) {
                    nextnum = 1;
                }
                respolygon.push_back(nextnum);
                regionsize++;
                if (nextnum == beginnum) {
                    break;
                }
                for (int i = 0; i < resize; i++) {
                    if (nextnum == repeatp1[i]) {
                        mark = false;
                        pointtmp1.pointx = outrings[region1].Ringx[nextnum];
                        pointtmp1.pointy = outrings[region1].Ringy[nextnum];
                        for (int j = 1; j < outrings[region2].Ringx.size(); j++) {
                            pointtmp2.pointx = outrings[region2].Ringx[j];
                            pointtmp2.pointy = outrings[region2].Ringy[j];
                            if (isequal(pointtmp1, pointtmp2)) {
                                nextnum = j;
                                break;
                            }
                        }
                        break;
                    }
                }
            }
            else {
                nextnum++;
                if (nextnum == outrings[region2].Ringx.size()) {
                    nextnum = 1;
                }
                respolygon.push_back(-nextnum);
                regionsize++;
                for (int i = 0; i < resize; i++) {
                    if (nextnum == repeatp2[i]) {
                        mark = true;
                        pointtmp1.pointx = outrings[region2].Ringx[nextnum];
                        pointtmp1.pointy = outrings[region2].Ringy[nextnum];
                        for (int j = 1; j < outrings[region1].Ringx.size(); j++) {
                            pointtmp2.pointx = outrings[region1].Ringx[j];
                            pointtmp2.pointy = outrings[region1].Ringy[j];
                            if (isequal(pointtmp1, pointtmp2)) {
                                nextnum = j;
                                break;
                            }
                        }
                        break;
                    }
                }
            }
        }
        cout << endl;
        for (int i = 0; i < regionsize; i++) {
            cout << respolygon[i] << "\t";
        }
    }
    else {
        cout << "the region is not contiguous" << endl;
        return 0;
    }
    //store into the respolygon 


    GDALDriver* pDriver = NULL;
    pDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");

    string shpfilename;

    //shpfilename = "data/" + itoa(region1) + itoa(region2) + "merge.shp"

    GDALDataset* pDatasetW = NULL;
    pDatasetW = pDriver->Create("data/merge.shp", 0, 0, 0, GDT_Unknown, NULL);

    OGRLayer* pLayerW = NULL;
    pLayerW = pDatasetW->CreateLayer("0", SpatialReference, GeometryType, NULL);

    OGRFieldDefn FieldTmp("id", OFTInteger);
    pLayerW->CreateField(&FieldTmp);

    pLayerW = pDatasetW->GetLayer(0);

    OGRFeature* pFeatureW = NULL;
    pFeatureW = OGRFeature::CreateFeature(pLayerW->GetLayerDefn());

    OGRPolygon PolygonTmp;
    OGRLinearRing RingTmp;

    for (int i = 0; i < iFeatureNum; i++) {
        if (i != region1 && i != region2) {
            OGRPolygon PolygonTmp;
            OGRLinearRing RingTmp;
            pFeatureW->SetField("id", i);
            for (int j = 0; j < outrings[i].Ringx.size(); j++) {
                RingTmp.addPoint(outrings[i].Ringx[j], outrings[i].Ringy[j]);
            }
            RingTmp.closeRings();
            PolygonTmp.addRing(&RingTmp);

            for (int j = 0; j < inholes[i].size(); j++) {
                RingTmp.empty();
                for (int k = 0; k < inholes[i][j].Ringx.size(); k++) {
                    RingTmp.addPoint(inholes[i][j].Ringx[k], inholes[i][j].Ringy[k]);
                }
                RingTmp.closeRings();
                PolygonTmp.addRing(&RingTmp);
            }
            pFeatureW->SetGeometry(&PolygonTmp);
            pLayerW->CreateFeature(pFeatureW);
        }
    }

    pFeatureW->SetField("id", iFeatureNum);
    iFeatureNum++;
    for (int i = 0; i < regionsize; i++) {
        if (respolygon[i] > 0) {
            RingTmp.addPoint(outrings[region1].Ringx[abs(respolygon[i])], outrings[region1].Ringy[abs(respolygon[i])]);
        }
        else {
            RingTmp.addPoint(outrings[region2].Ringx[abs(respolygon[i])], outrings[region2].Ringy[abs(respolygon[i])]);
        }
    }
    RingTmp.closeRings();
    PolygonTmp.addRing(&RingTmp);
    //write outrings
    for (int i = 0; i < inholes[region1].size(); i++) {
        RingTmp.empty();
        for (int j = 0; j < inholes[region1][i].Ringx.size(); j++) {
            RingTmp.addPoint(inholes[region1][i].Ringx[j], inholes[region1][i].Ringy[j]);
        }
        RingTmp.closeRings();
        PolygonTmp.addRing(&RingTmp);
    }
    //write inholes1
    for (int i = 0; i < inholes[region2].size(); i++) {
        RingTmp.empty();
        for (int j = 0; j < inholes[region2][i].Ringx.size(); j++) {
            RingTmp.addPoint(inholes[region2][i].Ringx[j], inholes[region2][i].Ringy[j]);
        }
        RingTmp.closeRings();
        PolygonTmp.addRing(&RingTmp);
    }
    //write inholes2
    pFeatureW->SetGeometry(&PolygonTmp);
    pLayerW->CreateFeature(pFeatureW);

    return 0;
}

bool isequal(Point a, Point b) {
    double dx = fabs(a.pointx - b.pointx);
    double dy = fabs(a.pointy - b.pointy);
    if (dx < 1e-2 && dy < 1e-2)
        return true;
    else
        return false;
}


// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
