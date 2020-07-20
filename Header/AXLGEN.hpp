#include "AGM.hpp"

int AXLGEN(const string &Geometry_filename, const string &AGL_output_file_TEMP, const string &AGM_output_file_TEMP,
           const double initialTime, const double terminalTime, const int timeStep, const double dt, int nIter) {

    ifstream GeoInfo;
    FILE *output = nullptr;
    Section *section = nullptr;
    Region *region = nullptr;
    Thing boundaryPTS; // intersect points
    Thing boundaryPTS_temp;
    Thing *boundaryPTS_region = nullptr;
    Thing_AxialLine XaxialLINE_temp;
    Thing_AxialLine YaxialLINE_temp;
    Thing *CrossPTS = nullptr;
    Thing *CrossPTS_temp = nullptr;
    Thing_AxialLine *XaxialLINE = nullptr;
    Thing_AxialLine *YaxialLINE = nullptr;
    int nD = 0; // dimension
    int nS = 0, nR = 0, nE = 0, nG = 0;// the number of sections, regions, elements, sign
    int PTSnum = 0; // Total number of points concluding cross points & boundary points
    int count = 0;
    string AGL_output_file;
    string AGM_output_file;
    stringstream ss;
    string num;

    for (size_t t = 0; t < nIter; t++) {

        PTSnum = 0;
        nS = 0, nR = 0, nE = 0, nG = 0;

        ss.str("");
        ss << t;
        num = ss.str();
        if (t > 1) {
            AGL_output_file = AGL_output_file_TEMP + num;
            AGM_output_file = AGM_output_file_TEMP + num;
        } else {
            AGL_output_file = AGL_output_file_TEMP;
            AGM_output_file = AGM_output_file_TEMP;
        }

        ReadInputFile(Geometry_filename, AGL_output_file, GeoInfo, output, &nD);

        switch (nD) {
            case D2:

                ReadGeoInfoNumber(GeoInfo, &nS, &nR, &nE, &nG);
                ReadGeoInfo(GeoInfo, section, region, nS, nR, nE, nD, nG, t);

                if (!CrossPTS) CrossPTS = new Thing[nR];
                if (!CrossPTS_temp) CrossPTS_temp = new Thing[nR];
                if (!boundaryPTS_region) boundaryPTS_region = new Thing[nR];
                if (!XaxialLINE) XaxialLINE = new Thing_AxialLine[nR];
                if (!YaxialLINE) YaxialLINE = new Thing_AxialLine[nR];

                boundaryPTS.TakeoffHeadCell();

                // FILE *fp; fp = fopen ("XaxialLine.dat","w");
                // FILE *fq; fq = fopen ("YaxialLine.dat","w");

                for (size_t k = 0; k < nR; k++) {

                    if (!strncmp("FIX", region[k].getName().c_str(), 3) && t > 0) {

                        printf("%s\n", "Find Fixed Region");

                        Write2DCrossPoints(output, CrossPTS[k], &region[k], &PTSnum);
                        Write2DBoundaryPoints(output, XaxialLINE[k], YaxialLINE[k], &region[k], &PTSnum,
                                              Count2DBoundaryPoints(output, XaxialLINE[k], YaxialLINE[k]));
                        Write2DXaxialline(output, CrossPTS[k], XaxialLINE[k], &region[k],
                                          PTSnum - 2 * XaxialLINE[k].Howmany() - 2 * YaxialLINE[k].Howmany());
                        Write2DYaxialline(output, CrossPTS[k], YaxialLINE[k], &region[k],
                                          PTSnum - 2 * YaxialLINE[k].Howmany());

                        boundaryPTS_region[k].GoHead();
                        boundaryPTS.GoTail();
                        boundaryPTS.Import(boundaryPTS_region[k], boundaryPTS_region[k].Howmany(), '+');

                        continue;
                    }
                    CrossPTS[k].TakeoffHeadCell();
                    CrossPTS_temp[k].TakeoffHeadCell();
                    XaxialLINE[k].TakeoffHeadLine();
                    YaxialLINE[k].TakeoffHeadLine();
                    count = 0;
                    printf("\n%s%s%s\n", "Region [", region[k].getName().c_str(), "] Start");
                    printf("\n%s\n", "Meeting background lines with boundary:");

                    for (size_t j = 0; j < region[k].getXgridnum(); j++) {
                        for (size_t i = 0; i < region[k].getSectionnum(); i++)
                            MasterMeetSlave(1, region[k].get2Dxgrid(j), region[k].getSection(i)->getEltnum(),
                                            region[k].getSection(i)->get2DElt(0), boundaryPTS_temp,
                                            region[k].getSign(i));
                        CopyTempToPerm(boundaryPTS, boundaryPTS_temp, XaxialLINE[k], XaxialLINE_temp, 'X');
                        if (!strncmp("FIX", region[k].getName().c_str(), 3)) {
                            CopyTempToPerm(boundaryPTS_region[k], boundaryPTS_temp, XaxialLINE[k], XaxialLINE_temp,
                                           'X');
                        }
                        printf("%s%d%s%d", "\r", ++count, "/", region[k].getXgridnum() + region[k].getYgridnum());
                    }

                    for (size_t j = 0; j < region[k].getYgridnum(); j++) {
                        for (size_t i = 0; i < region[k].getSectionnum(); i++)
                            MasterMeetSlave(1, region[k].get2Dygrid(j), region[k].getSection(i)->getEltnum(),
                                            region[k].getSection(i)->get2DElt(0), boundaryPTS_temp,
                                            region[k].getSign(i));
                        CopyTempToPerm(boundaryPTS, boundaryPTS_temp, YaxialLINE[k], YaxialLINE_temp, 'Y');
                        if (!strncmp("FIX", region[k].getName().c_str(), 3)) {
                            CopyTempToPerm(boundaryPTS_region[k], boundaryPTS_temp, YaxialLINE[k], YaxialLINE_temp,
                                           'Y');
                        }
                        printf("%s%d%s%d", "\r", ++count, "/", region[k].getXgridnum() + region[k].getYgridnum());
                    }
                    printf("\n%s\n", "All background lines met");

                    // XaxialLINE[k].GoHead();
                    // if (fp!=NULL) {
                    //   while (XaxialLINE[k].Run()) {
                    //     fprintf(fp, "%23.16e\t", XaxialLINE[k].Run()->Start()->Item(0));
                    //     fprintf(fp, "%23.16e\t", XaxialLINE[k].Run()->Start()->Item(1));
                    //     fprintf(fp, "%23.16e\t", XaxialLINE[k].Run()->Center()->getvalue(0));
                    //     fprintf(fp, "%23.16e\t", XaxialLINE[k].Run()->Center()->getvalue(1));
                    //     fprintf(fp, "%23.16e\t", XaxialLINE[k].Run()->End()->Item(0));
                    //     fprintf(fp, "%23.16e\n", XaxialLINE[k].Run()->End()->Item(1));
                    //     XaxialLINE[k].GoPost();
                    //   }
                    // }
                    //
                    // YaxialLINE[k].GoHead();
                    // if (fq!=NULL) {
                    //   while (YaxialLINE[k].Run()) {
                    //     fprintf(fq, "%23.16e\t", YaxialLINE[k].Run()->Start()->Item(0));
                    //     fprintf(fq, "%23.16e\t", YaxialLINE[k].Run()->Start()->Item(1));
                    //     fprintf(fq, "%23.16e\t", YaxialLINE[k].Run()->Center()->getvalue(0));
                    //     fprintf(fq, "%23.16e\t", YaxialLINE[k].Run()->Center()->getvalue(1));
                    //     fprintf(fq, "%23.16e\t", YaxialLINE[k].Run()->End()->Item(0));
                    //     fprintf(fq, "%23.16e\n", YaxialLINE[k].Run()->End()->Item(1));
                    //     YaxialLINE[k].GoPost();
                    //   }
                    // }

                    Delete2DAxialLine(&region[k], XaxialLINE[k]);
                    Delete2DAxialLine(&region[k], YaxialLINE[k]);

                    Make2DCrossPoints(&region[k], XaxialLINE[k], YaxialLINE[k], CrossPTS_temp[k], PTSnum);
                    CopyCrossPoints(CrossPTS[k], CrossPTS_temp[k], PTSnum);

                    Write2DCrossPoints(output, CrossPTS[k], &region[k], &PTSnum);
                    Write2DBoundaryPoints(output, XaxialLINE[k], YaxialLINE[k], &region[k], &PTSnum,
                                          Count2DBoundaryPoints(output, XaxialLINE[k], YaxialLINE[k], &region[k],
                                                                &PTSnum));
                    Write2DXaxialline(output, CrossPTS[k], XaxialLINE[k], &region[k],
                                      PTSnum - 2 * XaxialLINE[k].Howmany() - 2 * YaxialLINE[k].Howmany());
                    Write2DYaxialline(output, CrossPTS[k], YaxialLINE[k], &region[k],
                                      PTSnum - 2 * YaxialLINE[k].Howmany());
                }

                fclose(output);
                // fclose (fp);
                // fclose (fq);

                AGM(AGL_output_file, AGM_output_file, initialTime, terminalTime, timeStep, dt);

                break;
        }
    }

    return 0;
}
