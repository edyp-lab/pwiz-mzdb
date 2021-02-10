﻿/*
 * Original author: Nicholas Shulman <nicksh .at. u.washington.edu>,
 *                  MacCoss Lab, Department of Genome Sciences, UW
 *
 * Copyright 2009 University of Washington - Seattle, WA
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
using System.Collections.Generic;
using pwiz.Common.Collections;

namespace pwiz.Common.Chemistry
{
    public class IsotopeAbundances : ImmutableDictionary<string, MassDistribution>
    {
// ReSharper disable InconsistentNaming
        public static readonly IsotopeAbundances Default;
// ReSharper restore InconsistentNaming

        private IsotopeAbundances(IDictionary<string, MassDistribution> dictionary)
            : base(dictionary)
        {
        }

        public IsotopeAbundances SetAbundances(IDictionary<string, MassDistribution> abundances)
        {
            var dict = new Dictionary<string, MassDistribution>(this);
            foreach (var entry in abundances)
            {
                dict[entry.Key] = entry.Value;
            }
            return new IsotopeAbundances(dict);
        }

        public IsotopeAbundances SetAbundances(string element, MassDistribution massDistribution)
        {
            return SetAbundances(new Dictionary<string, MassDistribution> {{element, massDistribution}});
        }
        
        static IsotopeAbundances()
        {
            // These values may be easily updated to current NIST values using the tool at pwiz_tools\Skyline\Executables\ParseIsotopeAbundancesFromNIST

            // ReSharper disable LocalizableElement
            var defaults = new Dictionary<string, double[]>
            {
                // These values obtained 3/25/2019 from https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some
                {"H", new[]{1.00782503223,0.999885,2.01410177812,0.000115,}},
                {"He", new[]{3.0160293201,0.00000134,4.00260325413,0.99999866,}},
                {"Li", new[]{6.0151228874,0.0759,7.0160034366,0.9241,}},
                {"Be", new[]{9.012183065,1.0,}},
                {"B", new[]{10.01293695,0.199,11.00930536,0.801,}},
                {"C", new[]{12.0000000,0.9893,13.00335483507,0.0107,}},
                {"N", new[]{14.00307400443,0.99636,15.00010889888,0.00364,}},
                {"O", new[]{15.99491461957,0.99757,16.99913175650,0.00038,17.99915961286,0.00205,}},
                {"F", new[]{18.99840316273,1.0,}},
                {"Ne", new[]{19.9924401762,0.9048,20.993846685,0.0027,21.991385114,0.0925,}},
                {"Na", new[]{22.9897692820,1.0,}},
                {"Mg", new[]{23.985041697,0.7899,24.985836976,0.1000,25.982592968,0.1101,}},
                {"Al", new[]{26.98153853,1.0,}},
                {"Si", new[]{27.97692653465,0.92223,28.97649466490,0.04685,29.973770136,0.03092,}},
                {"P", new[]{30.97376199842,1.0,}},
                {"S", new[]{31.9720711744,0.9499,32.9714589098,0.0075,33.967867004,0.0425,35.96708071,0.0001,}},
                {"Cl", new[]{34.968852682,0.7576,36.965902602,0.2424,}},
                {"Ar", new[]{35.967545105,0.003336,37.96273211,0.000629,39.9623831237,0.996035,}},
                {"K", new[]{38.9637064864,0.932581,39.963998166,0.000117,40.9618252579,0.067302,}},
                {"Ca", new[]{39.962590863,0.96941,41.95861783,0.00647,42.95876644,0.00135,43.95548156,0.02086,45.9536890,0.00004,47.95252276,0.00187,}},
                {"Sc", new[]{44.95590828,1.0,}},
                {"Ti", new[]{45.95262772,0.0825,46.95175879,0.0744,47.94794198,0.7372,48.94786568,0.0541,49.94478689,0.0518,}},
                {"V", new[]{49.94715601,0.00250,50.94395704,0.99750,}},
                {"Cr", new[]{49.94604183,0.04345,51.94050623,0.83789,52.94064815,0.09501,53.93887916,0.02365,}},
                {"Mn", new[]{54.93804391,1.0,}},
                {"Fe", new[]{53.93960899,0.05845,55.93493633,0.91754,56.93539284,0.02119,57.93327443,0.00282,}},
                {"Co", new[]{58.93319429,1.0,}},
                {"Ni", new[]{57.93534241,0.68077,59.93078588,0.26223,60.93105557,0.011399,61.92834537,0.036346,63.92796682,0.009255,}},
                {"Cu", new[]{62.92959772,0.6915,64.92778970,0.3085,}},
                {"Zn", new[]{63.92914201,0.4917,65.92603381,0.2773,66.92712775,0.0404,67.92484455,0.1845,69.9253192,0.0061,}},
                {"Ga", new[]{68.9255735,0.60108,70.92470258,0.39892,}},
                {"Ge", new[]{69.92424875,0.2057,71.922075826,0.2745,72.923458956,0.0775,73.921177761,0.3650,75.921402726,0.0773,}},
                {"As", new[]{74.92159457,1.0,}},
                {"Se", new[]{73.922475934,0.0089,75.919213704,0.0937,76.919914154,0.0763,77.91730928,0.2377,79.9165218,0.4961,81.9166995,0.0873,}},
                {"Br", new[]{78.9183376,0.5069,80.9162897,0.4931,}},
                {"Kr", new[]{77.92036494,0.00355,79.91637808,0.02286,81.91348273,0.11593,82.91412716,0.11500,83.9114977282,0.56987,85.9106106269,0.17279,}},
                {"Rb", new[]{84.9117897379,0.7217,86.9091805310,0.2783,}},
                {"Sr", new[]{83.9134191,0.0056,85.9092606,0.0986,86.9088775,0.0700,87.9056125,0.8258,}},
                {"Y", new[]{88.9058403,1.0,}},
                {"Zr", new[]{89.9046977,0.5145,90.9056396,0.1122,91.9050347,0.1715,93.9063108,0.1738,95.9082714,0.0280,}},
                {"Nb", new[]{92.9063730,1.0,}},
                {"Mo", new[]{91.90680796,0.1453,93.90508490,0.0915,94.90583877,0.1584,95.90467612,0.1667,96.90601812,0.0960,97.90540482,0.2439,99.9074718,0.0982,}},
                {"Tc", new[]{98.0,1.0,}},
                {"Ru", new[]{95.90759025,0.0554,97.9052868,0.0187,98.9059341,0.1276,99.9042143,0.1260,100.9055769,0.1706,101.9043441,0.3155,103.9054275,0.1862,}},
                {"Rh", new[]{102.9054980,1.0,}},
                {"Pd", new[]{101.9056022,0.0102,103.9040305,0.1114,104.9050796,0.2233,105.9034804,0.2733,107.9038916,0.2646,109.90517220,0.1172,}},
                {"Ag", new[]{106.9050916,0.51839,108.9047553,0.48161,}},
                {"Cd", new[]{105.9064599,0.0125,107.9041834,0.0089,109.90300661,0.1249,110.90418287,0.1280,111.90276287,0.2413,112.90440813,0.1222,113.90336509,0.2873,115.90476315,0.0749,}},
                {"In", new[]{112.90406184,0.0429,114.903878776,0.9571,}},
                {"Sn", new[]{111.90482387,0.0097,113.9027827,0.0066,114.903344699,0.0034,115.90174280,0.1454,116.90295398,0.0768,117.90160657,0.2422,118.90331117,0.0859,119.90220163,0.3258,121.9034438,0.0463,123.9052766,0.0579,}},
                {"Sb", new[]{120.9038120,0.5721,122.9042132,0.4279,}},
                {"Te", new[]{119.9040593,0.0009,121.9030435,0.0255,122.9042698,0.0089,123.9028171,0.0474,124.9044299,0.0707,125.9033109,0.1884,127.90446128,0.3174,129.906222748,0.3408,}},
                {"I", new[]{126.9044719,1.0,}},
                {"Xe", new[]{123.9058920,0.000952,125.9042983,0.000890,127.9035310,0.019102,128.9047808611,0.264006,129.903509349,0.040710,130.90508406,0.212324,131.9041550856,0.269086,133.90539466,0.104357,135.907214484,0.088573,}},
                {"Cs", new[]{132.9054519610,1.0,}},
                {"Ba", new[]{129.9063207,0.00106,131.9050611,0.00101,133.90450818,0.02417,134.90568838,0.06592,135.90457573,0.07854,136.90582714,0.11232,137.90524700,0.71698,}},
                {"La", new[]{137.9071149,0.0008881,138.9063563,0.9991119,}},
                {"Ce", new[]{135.90712921,0.00185,137.905991,0.00251,139.9054431,0.88450,141.9092504,0.11114,}},
                {"Pr", new[]{140.9076576,1.0,}},
                {"Nd", new[]{141.9077290,0.27152,142.9098200,0.12174,143.9100930,0.23798,144.9125793,0.08293,145.9131226,0.17189,147.9168993,0.05756,149.9209022,0.05638,}},
                {"Pm", new[]{145.0,1.0,}},
                {"Sm", new[]{143.9120065,0.0307,146.9149044,0.1499,147.9148292,0.1124,148.9171921,0.1382,149.9172829,0.0738,151.9197397,0.2675,153.9222169,0.2275,}},
                {"Eu", new[]{150.9198578,0.4781,152.9212380,0.5219,}},
                {"Gd", new[]{151.9197995,0.0020,153.9208741,0.0218,154.9226305,0.1480,155.9221312,0.2047,156.9239686,0.1565,157.9241123,0.2484,159.9270624,0.2186,}},
                {"Tb", new[]{158.9253547,1.0,}},
                {"Dy", new[]{155.9242847,0.00056,157.9244159,0.00095,159.9252046,0.02329,160.9269405,0.18889,161.9268056,0.25475,162.9287383,0.24896,163.9291819,0.28260,}},
                {"Ho", new[]{164.9303288,1.0,}},
                {"Er", new[]{161.9287884,0.00139,163.9292088,0.01601,165.9302995,0.33503,166.9320546,0.22869,167.9323767,0.26978,169.9354702,0.14910,}},
                {"Tm", new[]{168.9342179,1.0,}},
                {"Yb", new[]{167.9338896,0.00123,169.9347664,0.02982,170.9363302,0.1409,171.9363859,0.2168,172.9382151,0.16103,173.9388664,0.32026,175.9425764,0.12996,}},
                {"Lu", new[]{174.9407752,0.97401,175.9426897,0.02599,}},
                {"Hf", new[]{173.9400461,0.0016,175.9414076,0.0526,176.9432277,0.1860,177.9437058,0.2728,178.9458232,0.1362,179.9465570,0.3508,}},
                {"Ta", new[]{179.9474648,0.0001201,180.9479958,0.9998799,}},
                {"W", new[]{179.9467108,0.0012,181.94820394,0.2650,182.95022275,0.1431,183.95093092,0.3064,185.9543628,0.2843,}},
                {"Re", new[]{184.9529545,0.3740,186.9557501,0.6260,}},
                {"Os", new[]{183.9524885,0.0002,185.9538350,0.0159,186.9557474,0.0196,187.9558352,0.1324,188.9581442,0.1615,189.9584437,0.2626,191.9614770,0.4078,}},
                {"Ir", new[]{190.9605893,0.373,192.9629216,0.627,}},
                {"Pt", new[]{189.9599297,0.00012,191.9610387,0.00782,193.9626809,0.3286,194.9647917,0.3378,195.96495209,0.2521,197.9678949,0.07356,}},
                {"Au", new[]{196.96656879,1.0,}},
                {"Hg", new[]{195.9658326,0.0015,197.96676860,0.0997,198.96828064,0.1687,199.96832659,0.2310,200.97030284,0.1318,201.97064340,0.2986,203.97349398,0.0687,}},
                {"Tl", new[]{202.9723446,0.2952,204.9744278,0.7048,}},
                {"Pb", new[]{203.9730440,0.014,205.9744657,0.241,206.9758973,0.221,207.9766525,0.524,}},
                {"Bi", new[]{208.9803991,1.0,}},
                {"Po", new[]{209.0,1.0,}},
                {"At", new[]{210.0,1.0,}},
                {"Rn", new[]{222.0,1.0,}},
                {"Fr", new[]{223.0,1.0,}},
                {"Ra", new[]{226.0,1.0,}},
                {"Ac", new[]{227.0,1.0,}},
                {"Th", new[]{232.0380558,1.0,}},
                {"Pa", new[]{231.0358842,1.0,}},
                {"U", new[]{234.0409523,0.000054,235.0439301,0.007204,238.0507884,0.992742,}},
                {"Np", new[]{237.0,1.0,}},
                {"Pu", new[]{244.0,1.0,}},
                {"Am", new[]{243.0613813,1.0,}},
                {"Cm", new[]{246.0672238,1.0,}},
                {"Bk", new[]{249.0749877,1.0,}},
                {"Cf", new[]{251.0795886,1.0,}},
                {"Es", new[]{252.082980,1.0,}},
                {"Fm", new[]{257.0951061,1.0,}},
                {"Md", new[]{260.10365,1.0,}},
                {"No", new[]{259.10103,1.0,}},
                {"Lr", new[]{262.10961,1.0,}},
                {"Rf", new[]{267.12179,1.0,}},
                {"Db", new[]{268.12567,1.0,}},
                {"Sg", new[]{271.13393,1.0,}},
                {"Bh", new[]{272.13826,1.0,}},
                {"Hs", new[]{270.13429,1.0,}},
                {"Mt", new[]{276.15159,1.0,}},
                {"Ds", new[]{281.16451,1.0,}},
                {"Rg", new[]{280.16514,1.0,}},
                {"Cn", new[]{285.17712,1.0,}},
                {"Nh", new[]{284.17873,1.0,}},
                {"Fl", new[]{289.19042,1.0,}},
                {"Mc", new[]{288.19274,1.0,}},
                {"Lv", new[]{293.20449,1.0,}},
                {"Ts", new[]{292.20746,1.0,}},
                {"Og", new[]{294.21392,1.0,}},

            };

            // These values have been used in Skyline for a very long time, changing them even a little causes test havoc
            var traditional_Skyline_CHONPS = new Dictionary<string, double[]>
            {
                {"H",new []{1.007825035,0.999855,2.014101779,0.000145,}},
                {"C",new []{12.0,0.98916,13.0033548378,0.01084,}},
                {"N",new []{14.003074,0.99633,15.0001088984,0.00366,}},
                {"O",new []{15.99491463,0.997576009706,16.9991315,0.000378998479,17.9991604,0.002044991815,}},
                {"P",new []{30.973762,1.0,}},
                {"S",new []{31.9720707,0.95021,32.971456,0.00745,33.967866,0.04221,35.96708,0.00013,}},
             };

            // ReSharper restore LocalizableElement

            var dict = new Dictionary<string, MassDistribution>();
            foreach (var entry in defaults)
            {
                var isotopes = new Dictionary<double, double>();
                var list = traditional_Skyline_CHONPS.ContainsKey(entry.Key)
                    ? traditional_Skyline_CHONPS[entry.Key]
                    : entry.Value;
                for (int i = 0; i < list.Length; i += 2)
                {
                    isotopes.Add(list[i], list[i + 1]);
                }
                dict.Add(entry.Key, MassDistribution.NewInstance(isotopes, 0, 0));
            }
            Default = new IsotopeAbundances(dict);
        }
    }
}