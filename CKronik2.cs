﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;

namespace BullseyeSharp
{
    public class CPep
    {
        public int charge;
        public float intensity;
        public double monoMass;
        public double basePeak;
        public double xCorr;
        public string mods;
    }

    public class CScan
    {
        public List<CPep> vPep;
        public int scanNum;
        public string file;
        public float rTime;

        public CScan()
        {
            scanNum = 0;
            file = "";
            rTime = 0;
            vPep = new List<CPep>();
        }
        public void clear() { vPep.Clear();}

    }

    public class CProfileData
    {
        public float intensity;
        public bool interpolated;
        public double monoMass;
        public float rTime;
        public int scanNum;
        public double xCorr;

        //Constructor
        public CProfileData()
        {
            interpolated = false;
            scanNum = 0;
            intensity = 0.0f;
            rTime = 0.0f;
            monoMass = 0.0;
            xCorr = 0.0;
        }
    }

    public class CPepProfile
    {
        public int lowScan;
        public int highScan;
        public int bestScan;
        public int charge;
        public int MS2Events;

        public float rTime;
        public float firstRTime;
        public float lastRTime;
        public float intensity;
        public float sumIntensity;

        public double monoMass;
        public double basePeak;
        public double xCorr;

        public string mods;
        public string gene;
        public string sequence;

        public List<CProfileData> profile;

        public CPepProfile()
        {
            profile = new List<CProfileData>();
        }

        public void setPoints(int i){
            profile.Clear();
            for (int j = 0; j < i; j++) profile.Add(new CProfileData());
        }
}

    public struct iTwo
    {
        public int scan;
        public int pep;
        public iTwo(int a, int b)
        {
            scan = a;
            pep = b;
        }
    }

    class CKronik2
    {
        private int iGapTol;
        private int iMatchTol;
        private double dPPMTol;
        private List<CPepProfile> vPeps;

        public CKronik2(){
            dPPMTol   = 10.0; 
            iGapTol   = 1;   
            iMatchTol = 3;
            vPeps = new List<CPepProfile>();
        }

        public CPepProfile this[int i]
        {
            get => vPeps[i];
            set => vPeps[i] = value;
        }

        public bool findMax(ref List<CScan> v, ref int s, ref int p)
        {
            float max = 0;
            int i;
            p = 0;
            for (i = 0; i < v.Count(); i++)
            {
                if (v[i].vPep.Count() == 0) continue;
                if (v[i].vPep[0].intensity > max)
                {
                    max = v[i].vPep[0].intensity;
                    s = i;
                }
            }
            return true;
        }

        public double interpolate(int x1, int x2, double y1, double y2, int x)
        {
            double m = (y2 - y1) / (x2 - x1);
            double b = y2 - m * x2;
            return m * x + b;
        }

        public bool processHK(string inFile, string outFile="")
        {
            char tag;
            int sIndex = 0;
            int pIndex = 0;
            int i, j, k, k1, k2;

            int pepCount = 0;
            List<CScan> allScans = new List<CScan>();

            double mass;
            double ppm;
            int charge;
            int gap;
            int matchCount;
            bool bMatch;

            //CPepProfile s;
            //CProfileData p;

            //for tracking which peptides
            List<iTwo> vLeft = new List<iTwo>();
            List<iTwo> vRight = new List<iTwo>();

            //clear data
            vPeps.Clear();

            //Read in the Hardklor results
            if (!System.IO.File.Exists(inFile)) return false;
            foreach (string line in System.IO.File.ReadLines(inFile))
            {
                tag = line[0];
                if (tag == 'S')
                {
                    allScans.Add(new CScan());
                    string[] tokens = line.Split('\t','\n','\r');
                    allScans.Last().scanNum = Convert.ToInt32(tokens[1]);
                    allScans.Last().rTime = (float)Convert.ToDouble(tokens[2]);
                    allScans.Last().file = tokens[3];
                }
                else
                {
                    pepCount++;
                    allScans.Last().vPep.Add(new CPep());
                    string[] tokens = line.Split('\t', '\n', '\r');
                    allScans.Last().vPep.Last().monoMass = Convert.ToDouble(tokens[1]);
                    allScans.Last().vPep.Last().charge = Convert.ToInt32(tokens[2]);
                    allScans.Last().vPep.Last().intensity = (float)Convert.ToDouble(tokens[3]);
                    allScans.Last().vPep.Last().basePeak = Convert.ToDouble(tokens[4]);
                    allScans.Last().vPep.Last().mods = tokens[7];
                    allScans.Last().vPep.Last().xCorr = Convert.ToDouble(tokens[8]);
                }
            }

            Console.WriteLine("" + pepCount + " peptides from " + allScans.Count() + " scans.\n");
            for (i = 0; i < allScans.Count; i++) allScans[i].vPep.Sort((a, b) => b.intensity.CompareTo(a.intensity));
            //for (i = 0; i < 10; i++) Console.WriteLine("" + i + "\t" + allScans[0].vPep[i].intensity);
            Console.WriteLine("Finding persistent peptide signals:");

            int startCount = pepCount;
            int lastPercent = 0;

            Console.Write(lastPercent);// cout << lastPercent;

            //Perform the Kronik analysis
            while (pepCount > 0)
            {
                if (!findMax(ref allScans, ref sIndex, ref pIndex)) break;

                mass = allScans[sIndex].vPep[pIndex].monoMass;
                charge = allScans[sIndex].vPep[pIndex].charge;
                matchCount = 1;

                //look left
                vLeft.Clear();
                gap = 0;
                i = sIndex - 1;
                while (i > -1 && gap <= iGapTol)
                {
                    bMatch = false;
                    int iScan = i;
                    int iPep = -1;
                    for (j = 0; j < allScans[i].vPep.Count; j++)
                    {
                        //if(allScans[i].vPep->at(j).intensity<0.0) continue;
                        ppm = (allScans[i].vPep[j].monoMass - mass) / mass * 1000000;
                        if (Math.Abs(ppm) < dPPMTol && allScans[i].vPep[j].charge == charge)
                        {
                            iPep = j;
                            gap = 0;
                            bMatch = true;
                            matchCount++;
                            break;
                        }
                    }
                    if (!bMatch) gap++;
                    vLeft.Add(new iTwo(iScan,iPep));
                    i--;
                }

                //look right
                vRight.Clear();
                gap = 0;
                i = sIndex + 1;
                while (i < allScans.Count && gap <= iGapTol)
                {
                    bMatch = false;
                    int iScan = i;
                    int iPep = -1;
                    for (j = 0; j < allScans[i].vPep.Count; j++)
                    {
                        //if(allScans[i].vPep->at(j).intensity<0.0) continue;
                        ppm = (allScans[i].vPep[j].monoMass - mass) / mass * 1000000;
                        if (Math.Abs(ppm) < dPPMTol && allScans[i].vPep[j].charge == charge)
                        {
                            iPep = j;
                            gap = 0;
                            bMatch = true;
                            matchCount++;
                            break;
                        }
                    }
                    if (!bMatch) gap++;
                    vRight.Add(new iTwo(iScan,iPep));
                    i++;
                }

                //Only keep sufficient matches
                if (matchCount >= iMatchTol)
                {
                    //trim gaps
                    while (vLeft.Count > 0 && vLeft[vLeft.Count - 1].pep < 0) vLeft.RemoveAt(vLeft.Count-1);
                    while (vRight.Count > 0 && vRight[vRight.Count - 1].pep < 0) vRight.RemoveAt(vRight.Count-1);

                    //apply basic information
                    vPeps.Add(new CPepProfile());
                    vPeps.Last().rTime = allScans[sIndex].rTime;
                    vPeps.Last().basePeak = allScans[sIndex].vPep[pIndex].basePeak;
                    vPeps.Last().bestScan = allScans[sIndex].scanNum;
                    vPeps.Last().charge = allScans[sIndex].vPep[pIndex].charge;
                    vPeps.Last().intensity = allScans[sIndex].vPep[pIndex].intensity;
                    vPeps.Last().monoMass = allScans[sIndex].vPep[pIndex].monoMass;
                    vPeps.Last().mods = allScans[sIndex].vPep[pIndex].mods;
                    vPeps.Last().xCorr = allScans[sIndex].vPep[pIndex].xCorr;
                    vPeps.Last().setPoints(vLeft.Count + vRight.Count + 1);
                    //vPeps.Last().profile = new List<CProfileData>();
                    if (vLeft.Count == 0)
                    {
                        vPeps.Last().lowScan = allScans[sIndex].scanNum;
                        vPeps.Last().firstRTime = allScans[sIndex].rTime;
                    }
                    else
                    {
                        vPeps.Last().lowScan = allScans[vLeft[vLeft.Count - 1].scan].scanNum;
                        vPeps.Last().firstRTime = allScans[vLeft[vLeft.Count - 1].scan].rTime;
                    }
                    if (vRight.Count == 0)
                    {
                        vPeps.Last().highScan = allScans[sIndex].scanNum;
                        vPeps.Last().lastRTime = allScans[sIndex].rTime;
                    }
                    else
                    {
                        vPeps.Last().highScan = allScans[vRight[vRight.Count - 1].scan].scanNum;
                        vPeps.Last().lastRTime = allScans[vRight[vRight.Count - 1].scan].rTime;
                    }

                    //apply datapoints
                    vPeps.Last().profile[0].intensity = allScans[sIndex].vPep[pIndex].intensity;
                    vPeps.Last().profile[0].interpolated = false;
                    vPeps.Last().profile[0].monoMass = allScans[sIndex].vPep[pIndex].monoMass;
                    vPeps.Last().profile[0].rTime = allScans[sIndex].rTime;
                    vPeps.Last().profile[0].scanNum = allScans[sIndex].scanNum;
                    vPeps.Last().profile[0].xCorr = allScans[sIndex].vPep[pIndex].xCorr;

                    i = 1;
                    j = 0;
                    while (j < vLeft.Count)
                    {
                        //Handle interpolated data
                        if (vLeft[j].pep < 0)
                        {
                            k = 0;
                            while (vLeft[j + k].pep < 0) k++;
                            k2 = j + k;
                            if (j == 0)
                            {
                                vPeps.Last().profile[i].intensity = (float)interpolate(vPeps.Last().profile[0].scanNum, allScans[vLeft[k2].scan].scanNum, (double)vPeps.Last().profile[0].intensity, (double)allScans[vLeft[k2].scan].vPep[vLeft[k2].pep].intensity, allScans[vLeft[j].scan].scanNum);
                                vPeps.Last().profile[i].monoMass = interpolate(vPeps.Last().profile[0].scanNum, allScans[vLeft[k2].scan].scanNum, vPeps.Last().profile[0].monoMass, allScans[vLeft[k2].scan].vPep[vLeft[k2].pep].monoMass, allScans[vLeft[j].scan].scanNum);
                            }
                            else
                            {
                                k = 0;
                                while (j + k >= 0 && vLeft[j + k].pep < 0) k--;
                                k1 = j + k;
                                if (k1 < 0)
                                {
                                    vPeps.Last().profile[i].intensity = (float)interpolate(vPeps.Last().profile[0].scanNum, allScans[vLeft[k2].scan].scanNum, (double)vPeps.Last().profile[0].intensity, (double)allScans[vLeft[k2].scan].vPep[vLeft[k2].pep].intensity, allScans[vLeft[j].scan].scanNum);
                                    vPeps.Last().profile[i].monoMass = interpolate(vPeps.Last().profile[0].scanNum, allScans[vLeft[k2].scan].scanNum, vPeps.Last().profile[0].monoMass, allScans[vLeft[k2].scan].vPep[vLeft[k2].pep].monoMass, allScans[vLeft[j].scan].scanNum);
                                }
                                else
                                {
                                    vPeps.Last().profile[i].intensity = (float)interpolate(allScans[vLeft[k1].scan].scanNum, allScans[vLeft[k2].scan].scanNum, (double)allScans[vLeft[k1].scan].vPep[vLeft[k1].pep].intensity, (double)allScans[vLeft[k2].scan].vPep[vLeft[k2].pep].intensity, allScans[vLeft[j].scan].scanNum);
                                    vPeps.Last().profile[i].monoMass = interpolate(allScans[vLeft[k1].scan].scanNum, allScans[vLeft[k2].scan].scanNum, allScans[vLeft[k1].scan].vPep[vLeft[k1].pep].monoMass, allScans[vLeft[k2].scan].vPep[vLeft[k2].pep].monoMass, allScans[vLeft[j].scan].scanNum);
                                }
                            }
                            vPeps.Last().profile[i].interpolated = true;
                            vPeps.Last().profile[i].rTime = allScans[vLeft[j].scan].rTime;
                            vPeps.Last().profile[i].scanNum = allScans[vLeft[j].scan].scanNum;
                            vPeps.Last().profile[i].xCorr = 0.0;
                        }
                        else
                        {
                            vPeps.Last().profile[i].intensity = allScans[vLeft[j].scan].vPep[vLeft[j].pep].intensity;
                            vPeps.Last().profile[i].interpolated = false;
                            vPeps.Last().profile[i].monoMass = allScans[vLeft[j].scan].vPep[vLeft[j].pep].monoMass;
                            vPeps.Last().profile[i].rTime = allScans[vLeft[j].scan].rTime;
                            vPeps.Last().profile[i].scanNum = allScans[vLeft[j].scan].scanNum;
                            vPeps.Last().profile[i].xCorr = allScans[vLeft[j].scan].vPep[vLeft[j].pep].xCorr;
                        }
                        i++;
                        j++;
                    }

                    j = 0;
                    while (j < vRight.Count)
                    {
                        //Handle interpolated data
                        if (vRight[j].pep < 0)
                        {
                            k = 0;
                            while (vRight[j + k].pep < 0) k++;
                            k2 = j + k;
                            if (j == 0)
                            {
                                vPeps.Last().profile[i].intensity = (float)interpolate(vPeps.Last().profile[0].scanNum, allScans[vRight[k2].scan].scanNum, vPeps.Last().profile[0].intensity, allScans[vRight[k2].scan].vPep[vRight[k2].pep].intensity, allScans[vRight[j].scan].scanNum);
                                vPeps.Last().profile[i].monoMass = interpolate(vPeps.Last().profile[0].scanNum, allScans[vRight[k2].scan].scanNum, vPeps.Last().profile[0].monoMass, allScans[vRight[k2].scan].vPep[vRight[k2].pep].monoMass, allScans[vRight[j].scan].scanNum);
                            }
                            else
                            {
                                k = 0;
                                while ((j + k > -1) && vRight[j + k].pep < 0) k--;
                                k1 = j + k;
                                if (k1 < 0)
                                {
                                    vPeps.Last().profile[i].intensity = (float)interpolate(vPeps.Last().profile[0].scanNum, allScans[vRight[k2].scan].scanNum, vPeps.Last().profile[0].intensity, allScans[vRight[k2].scan].vPep[vRight[k2].pep].intensity, allScans[vRight[j].scan].scanNum);
                                    vPeps.Last().profile[i].monoMass = interpolate(vPeps.Last().profile[0].scanNum, allScans[vRight[k2].scan].scanNum, vPeps.Last().profile[0].monoMass, allScans[vRight[k2].scan].vPep[vRight[k2].pep].monoMass, allScans[vRight[j].scan].scanNum);
                                }
                                else
                                {
                                    vPeps.Last().profile[i].intensity = (float)interpolate(allScans[vRight[k1].scan].scanNum, allScans[vRight[k2].scan].scanNum, allScans[vRight[k1].scan].vPep[vRight[k1].pep].intensity, allScans[vRight[k2].scan].vPep[vRight[k2].pep].intensity, allScans[vRight[j].scan].scanNum);
                                    vPeps.Last().profile[i].monoMass = interpolate(allScans[vRight[k1].scan].scanNum, allScans[vRight[k2].scan].scanNum, allScans[vRight[k1].scan].vPep[vRight[k1].pep].monoMass, allScans[vRight[k2].scan].vPep[vRight[k2].pep].monoMass, allScans[vRight[j].scan].scanNum);
                                }
                            }
                            vPeps.Last().profile[i].interpolated = true;
                            vPeps.Last().profile[i].rTime = allScans[vRight[j].scan].rTime;
                            vPeps.Last().profile[i].scanNum = allScans[vRight[j].scan].scanNum;
                            vPeps.Last().profile[i].xCorr = 0.0;
                        }
                        else
                        {
                            vPeps.Last().profile[i].intensity = allScans[vRight[j].scan].vPep[vRight[j].pep].intensity;
                            vPeps.Last().profile[i].interpolated = false;
                            vPeps.Last().profile[i].monoMass = allScans[vRight[j].scan].vPep[vRight[j].pep].monoMass;
                            vPeps.Last().profile[i].rTime = allScans[vRight[j].scan].rTime;
                            vPeps.Last().profile[i].scanNum = allScans[vRight[j].scan].scanNum;
                            vPeps.Last().profile[i].xCorr = allScans[vRight[j].scan].vPep[vRight[j].pep].xCorr;
                        }
                        i++;
                        j++;
                    }

                    //sort the profile
                    vPeps.Last().profile.Sort((a, b) => a.scanNum.CompareTo(b.scanNum));

                    //summed intensity (including interpolation)
                    vPeps.Last().sumIntensity = 0.0f;
                    for (i = 0; i < vPeps.Last().profile.Count; i++) vPeps.Last().sumIntensity += vPeps.Last().profile[i].intensity;

                    //Erase datapoints already used
                    for (i = 0; i < vLeft.Count; i++)
                    {
                        if (vLeft[i].pep < 0) continue;
                        allScans[vLeft[i].scan].vPep.RemoveAt(vLeft[i].pep);
                        pepCount--;
                    }
                    for (i = 0; i < vRight.Count; i++)
                    {
                        if (vRight[i].pep < 0) continue;
                        allScans[vRight[i].scan].vPep.RemoveAt(vRight[i].pep);
                        pepCount--;
                    }
                }

                //erase the one we're looking at
                allScans[sIndex].vPep.RemoveAt(pIndex);
                pepCount--;

                //update percent
                int iPercent = 100 - (int)((float)pepCount / (float)startCount * 100.0);
                if (iPercent > lastPercent)
                {
                    Console.Write("\b\b\b" + iPercent); // cout << "\b\b\b" << iPercent;
                    lastPercent = iPercent;
                }

            }
            Console.WriteLine("");

            if (outFile.Length>0){
                StreamWriter f = new StreamWriter(outFile);

                //Heading line
                f.Write("File\tFirst Scan\tLast Scan\tNum of Scans\tCharge\tMonoisotopic Mass\tBase Isotope Peak\t");
                f.Write("Best Intensity\tSummed Intensity\tFirst RTime\tLast RTime\tBest RTime\tBest Correlation\tModifications\n");

                for (i = 0; i < vPeps.Count; i++)
                {
                    f.Write("NULL\t");
                    f.Write(vPeps[i].lowScan + "\t");
                    f.Write(vPeps[i].highScan + "\t");
                    f.Write(vPeps[i].profile.Count + "\t");
                    f.Write(vPeps[i].charge + "\t");
                    f.Write(vPeps[i].monoMass + "\t");
                    f.Write(vPeps[i].basePeak + "\t");
                    f.Write(vPeps[i].intensity + "\t");
                    f.Write(vPeps[i].sumIntensity + "\t");
                    f.Write(vPeps[i].firstRTime + "\t");
                    f.Write(vPeps[i].lastRTime + "\t");
                    f.Write(vPeps[i].rTime + "\t");
                    f.Write(vPeps[i].xCorr + "\t");
                    f.WriteLine(vPeps[i].mods);
                }
                f.Close();
            }

            return true;
        }

        public void removeContam(double contam)
        {
            vPeps = vPeps.FindAll(a => (a.lastRTime - a.firstRTime) <= contam);
        }
        public void removeMass(double minMass, double maxMass)
        {
            vPeps = vPeps.FindAll(a => (a.monoMass >= minMass && a.monoMass <= maxMass));
        }

        public void setGapTol(int i) { iGapTol = i; }
        public void setMatchTol(int i) { iMatchTol = i; }
        public void setPPMTol(double d) { dPPMTol=d; }
        public int size() { return vPeps.Count; }
        public void sortBasePeak()
        {
            vPeps.Sort((a, b) => a.basePeak.CompareTo(b.basePeak));
        }
    }
}
