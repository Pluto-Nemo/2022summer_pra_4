// summer_pra_4_1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include <string>
#include "gdal201/gdal/include/ogrsf_frmts.h"
#include "gdal201/gdal/include/gdal_priv.h"
#pragma comment(lib, "gdal_i.lib")
using namespace std;

typedef struct OutRing {
    vector<double> Ringx;
    vector<double> Ringy;
}OutRing;
typedef struct DireLine {
    int startpointnum;
    int endpointnum;
}DireLine;
typedef struct Point {
    double pointx;
    double pointy;
    int linknum;
}Point;
typedef struct OutArc {
    vector<DireLine> lines;
    int startnodenum;
    int endnodenum;
    int size;
    int left;
    int right;
}OutArc;
typedef struct Region {
    vector<int> arcnum;
}Region;



bool isequal(Point a, Point b);

/*double caculateA(OutRing ring, vector<OutRing> inholes);
double caculateAS(OutRing ring);
double caculateP(OutRing ring);*/

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



    /// <summary>
    /// code for Q2_1
    Point pointtmp;
    vector<Point> points;
    int pointnum = 0;
    bool flag;
    for (int i = 0; i < outrings.size(); i++) {
        for (int j = 0; j < outrings[i].Ringx.size(); j++) {
            pointtmp.pointx = outrings[i].Ringx[j];
            pointtmp.pointy = outrings[i].Ringy[j];
            pointtmp.linknum = 0;
            flag = true;
            for (int k = 0; k < pointnum; k++) {
                if (isequal(pointtmp, points[k])) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                points.push_back(pointtmp);
                pointnum++;
                cout << points[pointnum - 1].pointx << "\t" << points[pointnum - 1].pointy << endl;
            }
        }
    }

    for (int i = 0; i < inholes.size(); i++) {
        for (int j = 0; j < inholes[i].size(); j++) {
            for (int k = 0; k < inholes[i][j].Ringx.size() - 1; k++) {
                pointtmp.pointx = inholes[i][j].Ringx[k];
                pointtmp.pointy = inholes[i][j].Ringy[k];
                pointtmp.linknum = 2;
                points.push_back(pointtmp);
                pointnum++;
                cout << points[pointnum - 1].pointx << "\t" << points[pointnum - 1].pointy << endl;
            }
        }
    }

    cout << pointnum << endl;
    cout << endl;

    DireLine sline;
    vector<DireLine> singlelines;
    int linenum = 0;

    for (int i = 0; i < outrings.size(); i++) {
        flag = true;
        for (int j = 0; j < outrings[i].Ringx.size() - 1; j++) {
            flag = true;
            for (int k = 0; k < pointnum; k++) {
                pointtmp.pointx = outrings[i].Ringx[j];
                pointtmp.pointy = outrings[i].Ringy[j];
                pointtmp.linknum = 0;
                if (isequal(pointtmp, points[k])) {
                    sline.startpointnum = k;
                    break;
                }
            }//to mark the num of sline's startpoint;
            for (int k = 0; k < pointnum; k++) {
                pointtmp.pointx = outrings[i].Ringx[j + 1];
                pointtmp.pointy = outrings[i].Ringy[j + 1];
                pointtmp.linknum = 0;
                if (isequal(pointtmp, points[k])) {
                    sline.endpointnum = k;
                    break;
                }
            }//to mark the num of sline's endpoint;
            for (int k = 0; k < linenum; k++) {
                if ((sline.endpointnum == singlelines[k].endpointnum && sline.startpointnum == singlelines[k].startpointnum) ||
                    (sline.endpointnum == singlelines[k].startpointnum && sline.startpointnum == singlelines[k].endpointnum)) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                singlelines.push_back(sline);
                cout << singlelines[linenum].startpointnum << "\t" << singlelines[linenum].endpointnum << endl;
                linenum++;
            }
        }
        for (int k = 0; k < pointnum; k++) {
            pointtmp.pointx = outrings[i].Ringx[outrings[i].Ringx.size() - 1];
            pointtmp.pointy = outrings[i].Ringy[outrings[i].Ringx.size() - 1];
            pointtmp.linknum = 0;
            if (isequal(pointtmp, points[k])) {
                sline.startpointnum = k;
                break;
            }
        }//to mark the num of sline's startpoint;
    }
    cout << linenum << endl;

    for (int i = 0; i < pointnum; i++) {
        points[i].linknum = 0;
    }
    for (int i = 0; i < linenum; i++) {
        points[singlelines[i].endpointnum].linknum++;
        points[singlelines[i].startpointnum].linknum++;
    }

    //struct all points and lines
    /// </summary>



    /// <summary>
    /// code for Q2_2

    //OutArc arctmp;
    vector<OutArc> outarcs;
    int outarcnum = 0;
    for (int i = 0; i < outrings.size(); i++) {
        int snode = 0;
        flag = true;

        for (int j = 0; j < outrings[i].Ringx.size(); j++) {
            pointtmp.pointx = outrings[i].Ringx[j];
            pointtmp.pointy = outrings[i].Ringy[j];
            for (int k = 0; k < pointnum; k++) {
                if (isequal(pointtmp, points[k]) && points[k].linknum > 2) {
                    flag = false;
                    snode = j;
                    break;
                }
            }
            if (!flag) {
                break;
            }
        }
        //find first node in i's outring whose num is snode


        int j = snode;
        flag = true;
        OutArc arctmp;
        while (j != snode || flag == true) {
            flag = false;
            int pn;
            pointtmp.pointx = outrings[i].Ringx[j];
            pointtmp.pointy = outrings[i].Ringy[j];
            for (int k = 0; k < pointnum; k++) {
                if (isequal(pointtmp, points[k])) {
                    pn = k;
                }
            }
            //find what the point number is (pn)

            bool reap = false;

            if (j == snode) {
                arctmp.startnodenum = pn;
                arctmp.lines.clear();
                arctmp.size = 0;
                sline.startpointnum = pn;
                arctmp.left = -1;
                arctmp.right = -1;
            }
            else if (points[pn].linknum > 2) {
                sline.endpointnum = pn;
                arctmp.lines.push_back(sline);
                arctmp.size++;
                sline.startpointnum = pn;
                arctmp.endnodenum = pn;

                //去重
                for (int k = 0; k < outarcnum; k++) {
                    reap = false;
                    if (outarcs[k].startnodenum == arctmp.startnodenum && outarcs[k].endnodenum == arctmp.endnodenum && outarcs[k].size == arctmp.size) {
                        reap = true;
                        outarcs[k].left = i;
                        break;
                    }
                    if (outarcs[k].startnodenum == arctmp.endnodenum && outarcs[k].endnodenum == arctmp.startnodenum && outarcs[k].size == arctmp.size) {
                        reap = true;
                        outarcs[k].left = i;
                        break;
                    }
                }
                if (reap == false) {
                    arctmp.right = i;
                    outarcs.push_back(arctmp);

                    cout << outarcnum << "\t" << outarcs[outarcnum].startnodenum << "\t" << outarcs[outarcnum].endnodenum << "\t"
                        << outarcs[outarcnum].right << "\t" << outarcs[outarcnum].left << endl;
                    for (int k = 0; k < outarcs[outarcnum].size; k++) {
                        cout << outarcs[outarcnum].lines[k].startpointnum << outarcs[outarcnum].lines[k].endpointnum << endl;
                    }
                    cout << endl;
                    //hint output

                    outarcnum++;
                }

                //new arc
                arctmp.startnodenum = pn;
                arctmp.lines.clear();
                arctmp.size = 0;
            }
            else {
                sline.endpointnum = pn;
                arctmp.lines.push_back(sline);
                arctmp.size++;
                sline.startpointnum = pn;
            }
            //push a new arc in outarcs;

            j++;
            if (j == outrings[i].Ringx.size() - 1) {
                j = 0;
            }
            //loop add j

            if (j == snode) {
                pointtmp.pointx = outrings[i].Ringx[j];
                pointtmp.pointy = outrings[i].Ringy[j];
                for (int k = 0; k < pointnum; k++) {
                    if (isequal(pointtmp, points[k])) {
                        pn = k;
                    }
                }
                sline.endpointnum = pn;
                arctmp.lines.push_back(sline);
                arctmp.size++;
                arctmp.endnodenum = pn;
                //去重
                for (int k = 0; k < outarcnum; k++) {
                    reap = false;
                    if (outarcs[k].startnodenum == arctmp.startnodenum && outarcs[k].endnodenum == arctmp.endnodenum && outarcs[k].size == arctmp.size) {
                        reap = true;
                        outarcs[k].left = i;
                        break;
                    }
                    if (outarcs[k].startnodenum == arctmp.endnodenum && outarcs[k].endnodenum == arctmp.startnodenum && outarcs[k].size == arctmp.size) {
                        reap = true;
                        outarcs[k].left = i;
                        break;
                    }
                }
                if (reap == false) {
                    arctmp.right = i;
                    outarcs.push_back(arctmp);

                    cout << outarcnum << "\t" << outarcs[outarcnum].startnodenum << "\t" << outarcs[outarcnum].endnodenum << "\t"
                        << outarcs[outarcnum].right << "\t" << outarcs[outarcnum].left << endl;
                    for (int k = 0; k < outarcs[outarcnum].size; k++) {
                        cout << outarcs[outarcnum].lines[k].startpointnum << outarcs[outarcnum].lines[k].endpointnum << endl;
                    }
                    cout << endl;
                    //hint output
                    outarcnum++;
                }
            }//deal with the last node
        }
        //count
    }

    for (int i = 0; i < inholes.size(); i++) {
        for (int j = 0; j < inholes[i].size(); j++) {
            OutArc arctmp;
            int pn = 0;
            pointtmp.pointx = inholes[i][j].Ringx[0];
            pointtmp.pointy = inholes[i][j].Ringy[0];
            for (int s = 0; s < pointnum; s++) {
                if (isequal(pointtmp, points[s])) {
                    pn = s;
                }
            }
            arctmp.startnodenum = pn;
            arctmp.lines.clear();
            arctmp.size = 0;
            sline.startpointnum = pn;
            arctmp.right = i;
            arctmp.left = -1;

            for (int k = 1; k < inholes[i][j].Ringx.size() - 1; k++) {
                pointtmp.pointx = inholes[i][j].Ringx[k];
                pointtmp.pointy = inholes[i][j].Ringy[k];
                for (int s = 0; s < pointnum; s++) {
                    if (isequal(pointtmp, points[s])) {
                        pn = s;
                    }
                }
                //get point's num
                sline.endpointnum = pn;
                arctmp.lines.push_back(sline);
                arctmp.size++;
                sline.startpointnum = pn;
            }

            pointtmp.pointx = inholes[i][j].Ringx[inholes[i][j].Ringx.size() - 1];
            pointtmp.pointy = inholes[i][j].Ringy[inholes[i][j].Ringx.size() - 1];
            for (int s = 0; s < pointnum; s++) {
                if (isequal(pointtmp, points[s])) {
                    pn = s;
                }
            }
            sline.endpointnum = pn;
            arctmp.lines.push_back(sline);
            arctmp.size++;
            sline.startpointnum = pn;
            arctmp.endnodenum = pn;
            outarcs.push_back(arctmp);

            cout << outarcnum << "\t" << outarcs[outarcnum].startnodenum << "\t" << outarcs[outarcnum].endnodenum << "\t"
                << outarcs[outarcnum].right << "\t" << outarcs[outarcnum].left << endl;
            for (int k = 0; k < outarcs[outarcnum].size; k++) {
                cout << outarcs[outarcnum].lines[k].startpointnum << outarcs[outarcnum].lines[k].endpointnum << endl;
            }

            cout << endl;
            outarcnum++;
        }
    }

    cout << outarcnum << endl;
    //out arc deal done!
    /// </summary>



    /// <summary>
    /// code for Q2_3
    Region regiontmp;
    vector<Region> regions;
    vector<vector<int>> arclist;
    int listlong[100];
    int numofregion = 0;
    for (int i = 0; i < 100; i++) {
        listlong[i] = 0;
    }
    for (int i = 0; i < outarcnum; i++) {
        if (outarcs[i].left != -1) {
            listlong[outarcs[i].left]++;
        }
        if (outarcs[i].right != -1) {
            listlong[outarcs[i].right]++;
        }
        if (numofregion < outarcs[i].left || numofregion < outarcs[i].right) {
            numofregion = max(outarcs[i].left, outarcs[i].right);
        }
    }
    for (int i = 0; i <= numofregion; i++) {
        vector<int> tmp;
        arclist.push_back(tmp);
    }
    for (int i = 0; i < outarcnum; i++) {
        if (outarcs[i].left != -1) {
            arclist[outarcs[i].left].push_back(i);
        }
        if (outarcs[i].right != -1) {
            arclist[outarcs[i].right].push_back(i);
        }
        if (numofregion < outarcs[i].left || numofregion < outarcs[i].right) {
            numofregion = max(outarcs[i].left, outarcs[i].right);
        }
    }

    bool regionlist[100][100];
    for (int i = 0; i <= numofregion; i++) {
        for (int j = 0; j < listlong[i]; j++) {
            regionlist[i][j] = true;
        }
    }

    for (int i = 0; i <= numofregion; i++) {
        //find first topo optimistic-direction arc
        int mark = 1;
        int beginnum = -1;
        int nextnum = -1;
        //int offset = 0;
        for (int j = 0; j < listlong[i]; j++) {
            if (outarcs[arclist[i][j]].right == i && outarcs[arclist[i][j]].endnodenum != outarcs[arclist[i][j]].startnodenum) {
                regiontmp.arcnum.push_back(arclist[i][j]);
                beginnum = outarcs[arclist[i][j]].startnodenum;
                nextnum = outarcs[arclist[i][j]].endnodenum;
                cout << i << "\tI\t" << arclist[i][j] << "\t";
                regionlist[i][j] = false;
                break;
            }
            /*else if (regionlist[i][j] && outarcs[arclist[i][j]].endnodenum == outarcs[arclist[i][j]].startnodenum) {
                regiontmp.arcnum.push_back(arclist[i][j]);
                regionlist[i][offset] = false; offset++;
                cout << i << "\tI\t" << arclist[i][j] << "\t";
            }*/
        }
        if (beginnum == -1) {
            for (int j = 0; j < listlong[i]; j++) {
                if (outarcs[arclist[i][j]].endnodenum != outarcs[arclist[i][j]].startnodenum) {
                    regiontmp.arcnum.push_back(-arclist[i][j]);
                    beginnum = outarcs[arclist[i][j]].startnodenum;
                    nextnum = outarcs[arclist[i][j]].endnodenum;
                    cout << -arclist[i][j] << "\t";
                    regionlist[i][j] = false;
                    mark = -1;
                    break;
                }
            }
        }
        //go through all the arc in the list and find the outring sequence
        while (true) {
            for (int j = 0; j < listlong[i]; j++) {
                if (regionlist[i][j]) {
                    if (outarcs[arclist[i][j]].startnodenum == nextnum) {
                        regiontmp.arcnum.push_back(arclist[i][j] * mark);
                        cout << arclist[i][j] * mark << "\t";
                        regionlist[i][j] = false;
                        nextnum = outarcs[arclist[i][j]].endnodenum;
                        break;
                    }
                    else if (outarcs[arclist[i][j]].endnodenum == nextnum) {
                        regiontmp.arcnum.push_back(-arclist[i][j] * mark);
                        cout << -arclist[i][j] * mark << "\t";
                        regionlist[i][j] = false;
                        nextnum = outarcs[arclist[i][j]].startnodenum;
                        break;
                    }
                }
            }
            if (nextnum == beginnum) {
                break;
            }
        }
        //pushback the outring sequence

        for (int j = 0; j < listlong[i]; j++) {
            if (outarcs[arclist[i][j]].endnodenum == outarcs[arclist[i][j]].startnodenum) {
                if (!regionlist[i][j]) {
                    regionlist[i][j] = true;
                }
                else {
                    regiontmp.arcnum.push_back(arclist[i][j]);
                    cout << arclist[i][j] << "\t";
                }
            }
        }
        cout << endl;
        regions.push_back(regiontmp);
        regiontmp.arcnum.clear();
    }

    cout << endl;

    ofstream outfile("Polygon_Arcs.txt", ios::trunc);
    outfile << "PolygonID:\tArcsID:" << endl;
    for (int i = 0; i <= numofregion; i++) {
        outfile << i << "\t|\t";
        for (int j = 0; j < listlong[i]; j++) {
            outfile << regions[i].arcnum[j] /* << " " << regionlist[i][j]*/ << "\t";//write to files, and for debugging
        }
        outfile << endl;
    }
    /// </summary>


    /*
    GDALDriver* pDriver = NULL;
    pDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");

    GDALDataset* pDatasetW = NULL;
    pDatasetW = pDriver->Create("data/node.shp", 0, 0, 0, GDT_Unknown, NULL);

    OGRLayer* pLayerW = NULL;
    pLayerW = pDatasetW->CreateLayer("0", SpatialReference, wkbPoint, NULL);//GeometryType is wkbpoint

    OGRFieldDefn FieldTmp("Degree", OFTInteger);
    pLayerW->CreateField(&FieldTmp);
    pLayerW = pDatasetW->GetLayer(0);

    OGRFeature* pFeatureW = NULL;
    pFeatureW = OGRFeature::CreateFeature(pLayerW->GetLayerDefn());

    OGRPoint PointTmp;

    for (int i = 0; i < pointnum; i++) {
        OGRPoint PointTmp;
        if (points[i].linknum > 2) {
            pFeatureW->SetField("Degree", points[i].linknum);
            PointTmp = OGRPoint(points[i].pointx, points[i].pointy);
            pFeatureW->SetGeometry(&PointTmp);
            pLayerW->CreateFeature(pFeatureW);
        }
    }*/
    //write node to shpfile <node>

    /*
    GDALDriver* pDriver = NULL;
    pDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");

    GDALDataset* pDatasetW = NULL;
    pDatasetW = pDriver->Create("data/arcs.shp", 0, 0, 0, GDT_Unknown, NULL);

    OGRLayer* pLayerW = NULL;
    pLayerW = pDatasetW->CreateLayer("0", SpatialReference, wkbLineString, NULL);//GeometryType is wkblinestring

    OGRFieldDefn FieldTmp1("Leftcode", OFTInteger);
    pLayerW->CreateField(&FieldTmp1);

    OGRFieldDefn FieldTmp2("Rightcode", OFTInteger);
    pLayerW->CreateField(&FieldTmp2);

    pLayerW = pDatasetW->GetLayer(0);

    OGRFeature* pFeatureW = NULL;
    pFeatureW = OGRFeature::CreateFeature(pLayerW->GetLayerDefn());

    int num = 0;

    for (int i = 0; i < outarcnum; i++) {
        OGRLineString LineStringTmp;
        pFeatureW->SetField("Leftcode", outarcs[i].left);
        pFeatureW->SetField("Rightcode", outarcs[i].right);
        for (int j = 0; j < outarcs[i].size; j++) {
            num = outarcs[i].lines[j].startpointnum;
            LineStringTmp.addPoint(points[num].pointx, points[num].pointy);
        }
        num = outarcs[i].lines[outarcs[i].size - 1].endpointnum;
        LineStringTmp.addPoint(points[num].pointx, points[num].pointy);
        pFeatureW->SetGeometry(&LineStringTmp);
        pLayerW->CreateFeature(pFeatureW);
    }
    */
    //write arcs to shpfile <arcs>

    /*
    GDALDriver* pDriver = NULL;
    pDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");

    GDALDataset* pDatasetW = NULL;
    pDatasetW = pDriver->Create("data/region.shp", 0, 0, 0, GDT_Unknown, NULL);

    OGRLayer* pLayerW = NULL;
    pLayerW = pDatasetW->CreateLayer("0", SpatialReference, GeometryType, NULL);//GeometryType is wkbpoint

    OGRFieldDefn FieldTmp("id", OFTInteger);
    pLayerW->CreateField(&FieldTmp);

    pLayerW = pDatasetW->GetLayer(0);

    OGRFeature* pFeatureW = NULL;
    pFeatureW = OGRFeature::CreateFeature(pLayerW->GetLayerDefn());


    for (int i = 0; i <= numofregion; i++) {
        OGRPolygon PolygonTmp;
        OGRLinearRing RingTmp;
        flag = false;
        pFeatureW->SetField("id", i);
        for (int j = 0; j < listlong[i]; j++) {
            int arclinesize = outarcs[abs(regions[i].arcnum[j])].size;
            if (!regionlist[i][j]) {
                flag = true;
                if (regions[i].arcnum[j] < 0) {
                    for (int k = 0; k < arclinesize; k++) {
                        RingTmp.addPoint(points[outarcs[abs(regions[i].arcnum[j])].lines[arclinesize - 1 - k].endpointnum].pointx,
                            points[outarcs[abs(regions[i].arcnum[j])].lines[arclinesize - 1 - k].endpointnum].pointy);
                    }
                   // RingTmp.addPoint(points[outarcs[abs(regions[i].arcnum[j])].lines[0].startpointnum].pointx,
                     //   points[outarcs[abs(regions[i].arcnum[j])].lines[0].startpointnum].pointy);
                }
                else {
                    for (int k = 0; k < arclinesize; k++) {
                        RingTmp.addPoint(points[outarcs[abs(regions[i].arcnum[j])].lines[k].startpointnum].pointx,
                            points[outarcs[abs(regions[i].arcnum[j])].lines[k].startpointnum].pointy);
                    }
                    //RingTmp.addPoint(points[outarcs[abs(regions[i].arcnum[j])].lines[arclinesize - 1].].pointx,
                      //  points[outarcs[abs(regions[i].arcnum[j])].lines[k].startpointnum].pointy);
                }
            }
        }
        if (flag) {
            RingTmp.closeRings();
            PolygonTmp.addRing(&RingTmp);
        }
        //outrings struct over

        for (int j = 0; j < listlong[i]; j++) {
            RingTmp.empty();
            if (regionlist[i][j]) {
                int arclinesize = outarcs[abs(regions[i].arcnum[j])].size;
                for (int k = 0; k < arclinesize; k++) {
                    RingTmp.addPoint(points[outarcs[abs(regions[i].arcnum[j])].lines[k].startpointnum].pointx,
                        points[outarcs[abs(regions[i].arcnum[j])].lines[k].startpointnum].pointy);
                }
                RingTmp.closeRings();
                PolygonTmp.addRing(&RingTmp);
            }
        }
        //inholes struct over
        pFeatureW->SetGeometry(&PolygonTmp);
        pLayerW->CreateFeature(pFeatureW);
    }*/
    //write arcs to shpfile <region>

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
