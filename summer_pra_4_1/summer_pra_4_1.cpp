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


double caculateA(OutRing ring, vector<OutRing> inholes);
double caculateAS(OutRing ring);
double caculateP(OutRing ring);

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
        iFieldType =  pField->GetType();
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

    GDALDriver* pDriver = NULL;
    pDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");

    GDALDataset* pDatasetW = NULL;
    pDatasetW = pDriver->Create("data/0.shp", 0, 0, 0, GDT_Unknown, NULL);

    OGRLayer* pLayerW = NULL;
    pLayerW = pDatasetW->CreateLayer("0", SpatialReference, GeometryType, NULL);

    for (int i = 0; i < iFieldNum; i++) {
        OGRFieldDefn FieldTmp(names[i].c_str(), types[i]);
        pLayerW->CreateField(&FieldTmp);
    }

    pLayerW = pDatasetW->GetLayer(0);

    OGRFeature* pFeatureW = NULL;
    pFeatureW = OGRFeature::CreateFeature(pLayerW->GetLayerDefn());


    for (int i = 0; i < iFeatureNum; i++) {
        OGRPolygon PolygonTmp;
        OGRLinearRing RingTmp;
        pFeatureW->SetField(names[0].c_str(), FieldsValues_list[i][0].c_str());
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
        pFeatureW->SetField(names[1].c_str(), caculateP(outrings[i]));
        pFeatureW->SetField(names[2].c_str(), caculateA(outrings[i], inholes[i]));
        //writein perimeter & area
        pFeatureW->SetGeometry(&PolygonTmp);
        pLayerW->CreateFeature(pFeatureW);
    }
    return 0;
}

double caculateP(OutRing ring) {
    double perimeter = 0, perimetertmp = 0;
    double dx, dy;
    for (int i = 0; i < ring.Ringx.size() - 1; i++) {
        dx = ring.Ringx[i + 1] - ring.Ringx[i];
        dy = ring.Ringy[i + 1] - ring.Ringy[i];
        perimetertmp = pow(dx * dx + dy * dy, 0.5);
        perimeter += perimetertmp;
    }
    dx = ring.Ringx[ring.Ringx.size() - 1] - ring.Ringx[0];
    dy = ring.Ringy[ring.Ringy.size() - 1] - ring.Ringy[0];
    perimetertmp = pow(dx * dx + dy * dy, 0.5);
    perimeter += perimetertmp;
    return perimeter;
}

double caculateAS(OutRing ring) {
    double area = 0, areatmp = 0;
    double a, b, h;
    for (int i = 0; i < ring.Ringx.size() - 1; i++) {
        a = ring.Ringy[i];
        b = ring.Ringy[i + 1];
        h = ring.Ringx[i + 1] - ring.Ringx[i];
        areatmp = (a + b) * h / 2;
        area += areatmp;
    }
    a = ring.Ringy[ring.Ringx.size() - 1];
    b = ring.Ringy[0];
    h = ring.Ringx[0] - ring.Ringx[ring.Ringx.size() - 1];
    areatmp = (a + b) * h / 2;
    area += areatmp;
    return area;
}

double caculateA(OutRing ring, vector<OutRing> inholes) {
    double area = 0, areatmp = 0;
    area = caculateAS(ring);
    for (int i = 0; i < inholes.size(); i++) {
        areatmp = caculateAS(inholes[i]);
        area -= areatmp;
    }
    return area;
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
