using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.SamplesConsole
{
    class RVEkanoninkhsGewmetriasBuilder
    {

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceKanonikhGewmetriaRveExampleParameters(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 95, //150, // diastaseis
                L02 = 95, //150,
                L03 = 95, //40,
                hexa1 = discr1 * subdiscr1,// diakritopoihsh
                hexa2 = discr1 * subdiscr1,
                hexa3 = discr3
            };

            grapheneSheetParameters gp;
            gp = new grapheneSheetParameters()
            {
                // parametroi shell
                E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                ni_shell = 0.0607, // stathera poisson
                elem1 = discr1_shell * subdiscr1_shell,
                elem2 = discr1_shell * subdiscr1_shell,
                L1 = 50,// nm  // DIORTHOSI 2 graphene sheets
                L2 = 50,// nm
                L3 = 112.5096153846, // nm
                a1_shell = 0, // nm
                tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                //parametroi cohesive epifaneias
                T_o_3 = 0.05,// Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.05,// Gpa
                D_o_1 = 0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
            return gpmp;
        }

        public static void HexaElementsOnlyRVEwithRenumbering(Model model, rveMatrixParameters mp, double[,] Dq, string renumberingVectorPath)
        {
            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Rve Matrix
            double E_disp = mp.E_disp; //Gpa
            double ni_disp = mp.ni_disp; // stather Poisson
            double L01 = mp.L01; // diastaseis
            double L02 = mp.L02;
            double L03 = mp.L03;
            int hexa1 = mp.hexa1;// diakritopoihsh
            int hexa2 = mp.hexa2;
            int hexa3 = mp.hexa3;
            // Perioxh parametroi Rve Matrix ews edw


            // Perioxh Gewmetria shmeiwn
            int nodeCounter = 0;

            int nodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            for (int h1 = 0; h1 < hexa1 + 1; h1++)
            {
                for (int h2 = 0; h2 < hexa2 + 1; h2++)
                {
                    for (int h3 = 0; h3 < hexa3 + 1; h3++)
                    {
                        nodeID = renumbering.GetNewNodeNumbering(RVEExamplesBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                        nodeCoordY = -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
                        nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                        nodeCounter++;
                    }
                }
            }
            // Perioxh Gewmetria shmeiwn ews edw

            //Perioxh Eisagwgh elements
            int elementCounter = 0;
            int subdomainID = 1;

            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = E_disp,
                PoissonRatio = ni_disp,
            };
            Element e1;
            int ElementID;
            int[] globalNodeIDforlocalNode_i = new int[8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based
                        globalNodeIDforlocalNode_i[0] = renumbering.GetNewNodeNumbering(RVEExamplesBuilder.Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[1] = renumbering.GetNewNodeNumbering(RVEExamplesBuilder.Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[2] = renumbering.GetNewNodeNumbering(RVEExamplesBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[3] = renumbering.GetNewNodeNumbering(RVEExamplesBuilder.Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[4] = renumbering.GetNewNodeNumbering(RVEExamplesBuilder.Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[5] = renumbering.GetNewNodeNumbering(RVEExamplesBuilder.Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[6] = renumbering.GetNewNodeNumbering(RVEExamplesBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[7] = renumbering.GetNewNodeNumbering(RVEExamplesBuilder.Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));

                        e1 = new Element()
                        {
                            ID = ElementID,
                            ElementType = new Hexa8NLRAM_1(material1, 3, 3, 3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                        };

                        for (int j = 0; j < 8; j++)
                        {
                            e1.NodesDictionary.Add(globalNodeIDforlocalNode_i[j], model.NodesDictionary[globalNodeIDforlocalNode_i[j]]);
                        }
                        model.ElementsDictionary.Add(e1.ID, e1);
                        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
                        elementCounter++;
                    }
                }
            }
            //Perioxh Eisagwgh elements

            //Tuple<int, int> nodeElementCounters = new Tuple<int, int>(nodeCounter, elementCounter);
            // change one tuple value
            //nodeElementCounters = new Tuple<int, int>(nodeElementCounters.Item1, elementCounter);
            // get one tuple value
            //elementCounter = nodeElementCounters.Item2;            
            //return nodeElementCounters;

            int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
            int f_komvoi_rve = kuvos;
            int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
            int komvos;
            //Dq = new double[9, 3 * p_komvoi_rve];
            for (int j = 0; j < p_komvoi_rve; j++)
            {
                komvos = renumbering.GetNewNodeNumbering(f_komvoi_rve + j + 1);
                Dq[0, 3 * j + 0] = model.NodesDictionary[komvos].X;
                Dq[1, 3 * j + 1] = model.NodesDictionary[komvos].Y;
                Dq[2, 3 * j + 2] = model.NodesDictionary[komvos].Z;
                Dq[3, 3 * j + 0] = model.NodesDictionary[komvos].Y;
                Dq[4, 3 * j + 1] = model.NodesDictionary[komvos].Z;
                Dq[5, 3 * j + 2] = model.NodesDictionary[komvos].X;
                Dq[6, 3 * j + 0] = model.NodesDictionary[komvos].Z;
                Dq[7, 3 * j + 1] = model.NodesDictionary[komvos].X;
                Dq[8, 3 * j + 2] = model.NodesDictionary[komvos].Y;
            }


        }

        public static void Reference2RVEExample10000withRenumberingwithInput(Model model)
        {
            double[,] Dq;
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop\Fxk_p_komvoi_rve.txt";
            string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop\o_xsunol_gs_{0}.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEkanoninkhsGewmetriasBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            

            int graphene_sheets_number = 10;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
            double[][] ekk_xyz = new double[graphene_sheets_number][];
            //double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1))- ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            RVEkanoninkhsGewmetriasBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);

            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample10000withRenumberingwithInput_2GrSh(Model model)
        {
            double[,] Dq;
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_2GrSh\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_2GrSh\Fxk_p_komvoi_rve.txt";
            string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_2GrSh\o_xsunol_gs_{0}.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEkanoninkhsGewmetriasBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;


            int graphene_sheets_number = 2;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
            double[][] ekk_xyz = new double[graphene_sheets_number][];
            //double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            RVEkanoninkhsGewmetriasBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);

            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample10000withRenumberingwithInput_1GrSh(Model model)
        {
            double[,] Dq;
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh\Fxk_p_komvoi_rve.txt";
            string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh\o_xsunol_gs_{0}.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEkanoninkhsGewmetriasBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;


            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
            double[][] ekk_xyz = new double[graphene_sheets_number][];
            //double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            RVEkanoninkhsGewmetriasBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);

            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample10000withRenumberingwithInput_1GrSh_sparse_v2(Model model)
        {
            double[,] Dq;
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh\Fxk_p_komvoi_rve.txt";
            string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh\o_xsunol_gs_{0}.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEkanoninkhsGewmetriasBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;


            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
            double[][] ekk_xyz = new double[graphene_sheets_number][];
            //double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            RVEkanoninkhsGewmetriasBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);

            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping_v2(model, hostGroup, embdeddedGroup);
        }

        public static void CreateRandomDataForGeom(int n_graphene_sheets, grapheneSheetParameters gp, rveMatrixParameters mp,double sigma_f)
        {
            //PROELEFSI REFERENCE_kanonikh_gewmetria\REF2_50_000_renu_new_multiple_algorithms_check_develop
            // 'ekk_xyz','phi1_ij','phi2_ij','rot_phi_1','rot_phi_2'

            // grapheneParameters
            int elem1 = gp.elem1;
            int elem2 = gp.elem2;
            double[] L1 = new double[n_graphene_sheets]; // = gp.L1;// nm
            double[] L2 = new double[n_graphene_sheets]; //gp.L2;// nm
            double[] L3 = new double[n_graphene_sheets]; //gp.L3; // nm
            for (int j = 0; j < n_graphene_sheets; j++)
            {
                L1[j] = gp.L1;//nm
                L2[j] = gp.L2; //nm
                L3[j] = gp.L3;
            }

                //matrixParameters
            double L01 = mp.L01; // diastaseis
            double L02 = mp.L02;
            double L03 = mp.L03;
            int hexa1 = mp.hexa1;// diakritopoihsh
            int hexa2 = mp.hexa2;
            int hexa3 = mp.hexa3;
            
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
            // Perioxh parametroi Rve Matrix ews edw


            double[][] ekk_xyz = new double [ n_graphene_sheets][];
            for (int j = 0; j < n_graphene_sheets; j++)
            {
                ekk_xyz[j] = new double[3];
                ekk_xyz[j][0] = -0.5 * L01 + 0.03 * L01 + ((L01 - 0.06 * L01) * rand());
                ekk_xyz[j][1] = -0.5 * L02 + 0.03 * L02 + ((L02 - 0.06 * L02) * rand());
                ekk_xyz[j][2] = -0.5 * L03 + 0.03 * L03 + ((L03 - 0.06 * L03) * rand());
            }

            double[] rot_phi_1 = new double[n_graphene_sheets];
            double[] rot_phi_2 = new double[n_graphene_sheets];
            for (int j = 0; j < n_graphene_sheets; j++)
            {
                rot_phi_1[j] = Math.PI * rand();
                rot_phi_2[j] = 0.5 * Math.PI * rand();
            }

            // create boxes
            double[][,] line_points = new double[n_graphene_sheets + 1][,]; //n_graphene_sheets +1 pou einai to RVE
            double[][,] line_segments = new double[n_graphene_sheets + 1][,];
            double[][,] pl_points = new double[n_graphene_sheets + 1][,];
            double[][,] vec1s = new double[n_graphene_sheets + 1][,];
            double[][,] vec2s = new double[n_graphene_sheets + 1][,];
            double[][,] perp_vec3s = new double[n_graphene_sheets + 1][,];
            for (int j = 0; j < n_graphene_sheets+1; j++)
            {
                line_points[j] = new double[3, 12];
                line_segments[j] = new double[3, 12];
                pl_points[j] = new double[3, 6];
                vec1s[j] = new double[3, 6];
                vec2s[j] = new double[3, 6];
                perp_vec3s[j] = new double[3, 6];
            }

            //prwta ta n_graphene_sheets
            for (int j = 0; j < n_graphene_sheets ; j++)
            {
                createGrShBoxSurroundingPlanes(L1[j], L2[j], sigma_f, rot_phi_1[j], rot_phi_2[j],ekk_xyz[j], line_points[j], line_segments[j], pl_points[j], vec1s[j], vec2s[j], perp_vec3s[j]);
            }
            //kai meta to RVE box 
            createGrShBoxSurroundingPlanes(L01, L02, L03 / 7, 0,0,new double [3] { 0,0,0}, line_points[n_graphene_sheets], line_segments[n_graphene_sheets], pl_points[n_graphene_sheets], vec1s[n_graphene_sheets], vec2s[n_graphene_sheets], perp_vec3s[n_graphene_sheets]);



            bool[,] GrSh_intersections = new bool [n_graphene_sheets + 1, n_graphene_sheets + 1];
            for (int i1 = 0; i1 < n_graphene_sheets+1; i1++)
            { for (int j1 = 0; j1 < n_graphene_sheets+1; j1++) { GrSh_intersections[i1, j1] = false; } }  // arxiko ola einai false

            for (int i1 = 0; i1 < n_graphene_sheets + 1; i1++)
            {
                for (int i2 = i1+1; i2 < n_graphene_sheets + 1; i2++)
                {
                    //check for intersection between i1 and i2 graphene sheets
                    GrSh_intersections[i1, i2] = check_for_intersection_between_i1_and_i2_graphene_sheets(i1, i2, line_points, line_segments, pl_points, vec1s, vec2s, perp_vec3s);
                }
            }
            //Gemisma summetrikou
            for (int i1 = 0; i1 < n_graphene_sheets + 1; i1++)
            {
                for (int i2 = i1 + 1; i2 < n_graphene_sheets + 1; i2++)
                {
                    GrSh_intersections[i2,i1] = GrSh_intersections[i1, i2];
                }
            }

            // TODO display rve geometry
            //

            while (combine_bool_any(GrSh_intersections))
            {

                //epilogh GrSh
                double inter_counter_min = sum_of_bool_array_elements(GrSh_intersections);

                int chosen_GrSh=-1;
                for (int i1 = 0; i1 < n_graphene_sheets ; i1++)
                {
                    bool[,] new_intersections = Clone_array(GrSh_intersections);

                    for (int j = 0; j < new_intersections.GetLength(1); j++)
                    { new_intersections[i1, j] = false; }
                    for (int j = 0; j < new_intersections.GetLength(0); j++)
                    { new_intersections[j,i1] = false; }

                    double inter_counter = sum_of_bool_array_elements(new_intersections);

                    if (inter_counter<inter_counter_min)
                    {
                        chosen_GrSh = i1;
                        inter_counter_min = inter_counter;
                    }
                }


                // "afairesh" tou epilegmenou kai topothetisi neou
                rot_phi_1[chosen_GrSh] = Math.PI * rand();
                rot_phi_2[chosen_GrSh] =0.5* Math.PI * rand();
                ekk_xyz[chosen_GrSh][0] = -0.5 * L01 + 0.03 * L01 + ((L01 - 0.06 * L01) * rand());
                ekk_xyz[chosen_GrSh][1] = -0.5 * L02 + 0.03 * L02 + ((L02 - 0.06 * L02) * rand());
                ekk_xyz[chosen_GrSh][2] = -0.5 * L03 + 0.03 * L03 + ((L03 - 0.06 * L03) * rand());
                createGrShBoxSurroundingPlanes(L1[chosen_GrSh], L2[chosen_GrSh], sigma_f, rot_phi_1[chosen_GrSh], rot_phi_2[chosen_GrSh], ekk_xyz[chosen_GrSh], line_points[chosen_GrSh], line_segments[chosen_GrSh], pl_points[chosen_GrSh], vec1s[chosen_GrSh], vec2s[chosen_GrSh], perp_vec3s[chosen_GrSh]);


                for (int i1 = 0; i1 < n_graphene_sheets+1; i1++)
                {
                    //check for intersection between i1 and chosen_GrSh graphene sheets
                    if (i1==chosen_GrSh)
                    { }
                    else
                    {
                        GrSh_intersections[i1, chosen_GrSh] = check_for_intersection_between_i1_and_i2_graphene_sheets(i1, chosen_GrSh, line_points, line_segments, pl_points, vec1s, vec2s, perp_vec3s);
                    }
                }

                //Gemisma summetrikou
                for (int i1 = 0; i1 < n_graphene_sheets + 1; i1++)
                {
                    GrSh_intersections[ chosen_GrSh,i1] = GrSh_intersections[i1, chosen_GrSh];
                }
            }

            // return 

        }

        public static void createGrShBoxSurroundingPlanes(double L1,double L2,double sigma_f,double rot_phi_1,double rot_phi_2,double [] ekk_xyz,double[,] line_points, double[,] line_segments, double[,] pl_points, double[,] vec1s, double[,] vec2s, double[,]perp_vec3s)
        {
            double[,] line_points_data = new double[3, 12] { { 0.5 * L1, -0.5 * L1, -0.5 * L1, 0.5 * L1, -0.5 * L1, -0.5 * L1, -0.5 * L1, -0.5 * L1, 0.5 * L1, -0.5 * L1, -0.5 * L1, 0.5 * L1 },
                {0.5*L2,0.5*L2,-0.5*L2,-0.5*L2,0.5*L2,-0.5*L2,-0.5*L2,+0.5*L2, -0.5*L2,-0.5*L2,-0.5*L2,-0.5*L2 },
                { -3.5*sigma_f,-3.5*sigma_f,-3.5*sigma_f,-3.5*sigma_f,+3.5*sigma_f,+3.5*sigma_f,-3.5*sigma_f,-3.5*sigma_f,+3.5*sigma_f,+3.5*sigma_f,-3.5*sigma_f,-3.5*sigma_f} };
            double[,] line_segments_data = new double[3, 12] { {0,0,0,0,L1,L1,L1,L1,0,0,0,0 },
                {0,0,0,0,0,0,0,0,L2,L2,L2,L2 },
                {2*3.5*sigma_f,2*3.5*sigma_f,2*3.5*sigma_f,2*3.5*sigma_f,0,0,0,0,0,0,0,0 } };
            double[,] pl_points_data = new double[3, 6] { { -0.5*L1,-0.5*L1,-0.5*L1,+0.5*L1,+0.5*L1,+0.5*L1},
            { -0.5*L2,-0.5*L2,-0.5*L2,0.5*L2,0.5*L2,0.5*L2},
            { -3.5*sigma_f,-3.5*sigma_f,-3.5*sigma_f,3.5*sigma_f,3.5*sigma_f,3.5*sigma_f,} };
            double[,] vec1s_data = new double[3, 6] { { L1,0,0,-L1,0,0},
            {0,L2,0,0,-L2,0 },
            {0,0,2*3.5*sigma_f,0,0,-2*3.5*sigma_f }};
            double[,] vec2s_data = new double[3, 6] { {0,0,L1,0,0,-L1},
            {L2,0,0,-L2,0,0 },
            { 0,2*3.5*sigma_f,0,0,-2*3.5*sigma_f,0}};
            double[,] perp_vec3s_data = new double[3, 6] { {0,L1,0,0,-L1,0 },
            { 0,0,L2,0,0,-L2},
            { 2*3.5*sigma_f,0,0,-2*3.5*sigma_f,0,0}};
            for ( int i1=0;i1<3; i1++)
            {
                for (int j = 0; j <12; j++)
                {
                    line_points[i1, j] = line_points_data[i1, j];
                    line_segments[i1, j] = line_segments_data[i1, j];
                }
            }
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j = 0; j < 6; j++)
                {
                    pl_points[i1, j] = pl_points_data[i1, j];
                    vec1s[i1, j] = vec1s_data[i1, j];
                    vec2s[i1, j] = vec2s_data[i1, j];
                    perp_vec3s[i1, j] = perp_vec3s_data[i1, j];
                }
            }


            // ROTATION copy apo to gewmetria_shell_coh...._ekk_random th dhmiourgia tou
            // Qij kai meta efarmogh sta dianusmata(theshs kai katefthunshs)
            double e1_new_z = Math.Sin(rot_phi_2);
            double e1_new_y = Math.Sin(rot_phi_1) * Math.Cos(rot_phi_2); // e1_new_xy = cos(rot_phi_2);
            double e1_new_x = Math.Cos(rot_phi_1) * Math.Cos(rot_phi_2);

            double e2_new_y = Math.Cos(rot_phi_1);
            double e2_new_x = -Math.Sin(rot_phi_1);
            double e2_new_z = 0;

            double[,] e_new = new double[3, 3] { {e1_new_x,e2_new_x,0 }, {e1_new_y,e2_new_y,0 }, {e1_new_z,e2_new_z,0 } };
            double[] e_new_cross = new double[3];
            RVEExamplesBuilder.cross(new double[3] { e1_new_x, e1_new_y, e1_new_z }, new double[3] { e2_new_x, e2_new_y, e2_new_z }, e_new_cross);
            for (int i1 = 0; i1 < 3; i1++)
            { e_new[i1, 2] = e_new_cross[i1]; }

            double[,] e_old = new double[3, 3] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

            double[,] Qij = new double[3, 3];
            for (int q1 = 0; q1 < 3; q1++)
            {
                for (int q2 = 0; q2 < 3; q2++)
                {
                    Qij[q1, q2] = RVEExamplesBuilder.dot_product(new double[3] { e_old[0, q1], e_old[1, q1], e_old[2, q1] }, new double[3] { e_new[0, q2], e_new[1, q2], e_new[2, q2] });
                }
            }

            for (int q1 = 0; q1 < 12; q1++)
            {
                double[] product;
                product = RVEExamplesBuilder.MatVecMult(Qij, new double[3] { line_points[0, q1], line_points[1, q1], line_points[2, q1] });
                for (int q2 = 0; q2 < 3; q2++) { line_points[q2, q1] = product[q2]; }
                product = RVEExamplesBuilder.MatVecMult(Qij, new double[3] { line_segments[0, q1], line_segments[1, q1], line_segments[2, q1] });
                for (int q2 = 0; q2 < 3; q2++) { line_segments[q2, q1] = product[q2]; }
            }

            for (int q1 = 0; q1 < 6; q1++)
            {
                double[] product;
                product = RVEExamplesBuilder.MatVecMult(Qij, new double[3] { pl_points[0, q1], pl_points[1, q1], pl_points[2, q1] });
                for (int q2 = 0; q2 < 3; q2++) { pl_points[q2, q1] = product[q2]; }
                product = RVEExamplesBuilder.MatVecMult(Qij, new double[3] { vec1s[0, q1], vec1s[1, q1], vec1s[2, q1] });
                for (int q2 = 0; q2 < 3; q2++) { vec1s[q2, q1] = product[q2]; }
                product = RVEExamplesBuilder.MatVecMult(Qij, new double[3] { vec2s[0, q1], vec2s[1, q1], vec2s[2, q1] });
                for (int q2 = 0; q2 < 3; q2++) { vec2s[q2, q1] = product[q2]; }
                product = RVEExamplesBuilder.MatVecMult(Qij, new double[3] { perp_vec3s[0, q1], perp_vec3s[1, q1], perp_vec3s[2, q1] });
                for (int q2 = 0; q2 < 3; q2++) { perp_vec3s[q2, q1] = product[q2]; }
            }


            // TRANSLATION copy apo to gewmetria_shell_coh...._ekk_random kai to
            //efarmozoume mono sta points kai oxi kai sta dianusmata katefthunshs
            for (int q1 = 0; q1 < 12; q1++)
            {
                line_points[0, q1] += ekk_xyz[0];
                line_points[1, q1] += ekk_xyz[1];
                line_points[2, q1] += ekk_xyz[2];
            }
            for (int q1 = 0; q1 < 6; q1++)
            {
                pl_points[0, q1] += ekk_xyz[0];
                pl_points[1, q1] += ekk_xyz[1];
                pl_points[2, q1] += ekk_xyz[2];
            }
        }

        public static bool check_for_intersection_between_i1_and_i2_graphene_sheets(int i1, int i2, double[][,] line_points, double[][,] line_segments, double[][,] pl_points, double[][,] vec1s, double[][,] vec2s, double[][,] perp_vec3s)
        {
            bool they_intersect = false;

            bool stop_q1 = false;
            for (int q1 = 0; q1 < line_points[0].GetLength(1); q1++)
            {
                for (int q2 = 0; q2 < pl_points[0].GetLength(1); q2++)
                {
                    if (check_if_plane_and_vec_intersect(new double[3] { pl_points[i1][0, q2], pl_points[i1][1, q2], pl_points[i1][2, q2] },
                        new double[3] { vec1s[i1][0, q2], vec1s[i1][1, q2], vec1s[i1][2, q2] },
                        new double[3] { vec2s[i1][0, q2], vec2s[i1][1, q2], vec2s[i1][2, q2] },
                        new double[3] { perp_vec3s[i1][0, q2], perp_vec3s[i1][1, q2], perp_vec3s[i1][2, q2] },
                        new double[3] { line_points[i2][0, q1], line_points[i2][1, q1], line_points[i2][2, q1] },
                        new double[3] { line_segments[i2][0, q1], line_segments[i2][1, q1], line_segments[i2][2, q1] }))
                    {
                        they_intersect = true;
                        stop_q1 = true;
                        break;
                    }
                }
                if (stop_q1)
                {
                    break;
                }
            }

            if (stop_q1)
            {
            }
            else
            {
                for (int q1 = 0; q1 < line_points[0].GetLength(1); q1++)
                {
                    for (int q2 = 0; q2 < pl_points[0].GetLength(1); q2++)
                    {
                        if (check_if_plane_and_vec_intersect(new double[3] { pl_points[i2][0, q2], pl_points[i2][1, q2], pl_points[i2][2, q2] },
                            new double[3] { vec1s[i2][0, q2], vec1s[i2][1, q2], vec1s[i2][2, q2] },
                            new double[3] { vec2s[i2][0, q2], vec2s[i2][1, q2], vec2s[i2][2, q2] },
                            new double[3] { perp_vec3s[i2][0, q2], perp_vec3s[i2][1, q2], perp_vec3s[i2][2, q2] },
                            //
                            new double[3] { line_points[i1][0, q1], line_points[i1][1, q1], line_points[i1][2, q1] },
                            new double[3] { line_segments[i1][0, q1], line_segments[i1][1, q1], line_segments[i1][2, q1] }))
                        {
                            they_intersect = true;
                            stop_q1 = true;
                            break;
                        }
                    }
                    if (stop_q1)
                    {
                        break;
                    }
                }
            }

            return they_intersect;
        }

        public static bool check_if_plane_and_vec_intersect(double[] pl_point, double[] vec1, double[] vec2, double[] perp_vec3, double[] line_point, double[] line_segment )
        {
            bool intersection = false;

            bool plane_and_line_parallel = false;

            if ( RVEExamplesBuilder.dot_product(perp_vec3,line_segment)==0)
            {
                plane_and_line_parallel = true;
            }
            else
            {
                // diadikasia evreshs shmeiou tomhs
                double t = (perp_vec3[1] * (pl_point[1] - line_point[1]) + perp_vec3[2] * (pl_point[2] - line_point[2]) + perp_vec3[3] * (pl_point[3] - line_point[3])) /(perp_vec3[1] * line_segment[1] + perp_vec3[2] * line_segment[2] + perp_vec3[3] * line_segment[3]);
                double[] intersection_point = new double[3];
                for (int q2 = 0; q2 < 3; q2++)
                { intersection_point[q2] = line_point[q2] +t*line_segment[q2]; }
                double[] pl_point_TO_int_point_VECTOR = new double[3];
                for (int q2 = 0; q2 < 3; q2++)
                { pl_point_TO_int_point_VECTOR[q2] = intersection_point[q2] -pl_point[q2]; }

                double[,] A = new double[3, 3] { {vec1[0],vec2[0],perp_vec3[0] }, { vec1[1], vec2[1], perp_vec3[1] }, { vec1[2], vec2[2], perp_vec3[2] } };
                double[] lamda = invert3by3(A, pl_point_TO_int_point_VECTOR);

                if (t<=1)
                {
                    if (t>=0)
                    {
                        if (lamda[0] <=1)
                        {
                            if (lamda[0] >=0)
                            {
                                if (lamda[1] <= 1)
                                {
                                    if (lamda[1] >= 0)
                                    {
                                        intersection = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //PART B
            if(plane_and_line_parallel)
            {
                double distance = point_plane_distance(pl_point, perp_vec3, line_point);

                if (distance==0)
                {
                    bool[] plane_section = new bool[4];
                    plane_section[0] = check_if_vec_and_vec_intersect(pl_point, vec1, line_point, line_segment);
                    plane_section[1] = check_if_vec_and_vec_intersect(pl_point, vec2, line_point, line_segment);
                    plane_section[2] = check_if_vec_and_vec_intersect(new double[3] { pl_point[0] + vec2[0], pl_point[1] + vec2[1], pl_point[2] + vec2[2] }, vec1, line_point, line_segment);
                    plane_section[3] = check_if_vec_and_vec_intersect(new double[3] { pl_point[0] + vec1[0],pl_point[1] + vec1[1], pl_point[2] + vec1[2] }, vec2, line_point, line_segment);

                    if(plane_section[0]|| plane_section[0]|| plane_section[0]|| plane_section[0])
                    {
                        intersection = true;
                    }

                }                
            }
            return intersection;
        }

        public static double point_plane_distance(double [] pl_point, double[] perp_vec3, double[] line_point)
        {
            double[] QP = new double[3];
            for (int q2 = 0; q2 < 3; q2++)
            { QP[q2] = line_point[q2] - pl_point[q2]; }

            double[] n = new double[3];
            double norm_perp_vec3 = Math.Sqrt(perp_vec3[0] * perp_vec3[0] + perp_vec3[1] * perp_vec3[1] + perp_vec3[2] * perp_vec3[2]);
            for (int q2 = 0; q2 < 3; q2++)
            { n[q2] = perp_vec3[q2] / norm_perp_vec3; }

            double distance = RVEExamplesBuilder.dot_product(QP, n);

            return distance;
        }

        public static bool check_if_vec_and_vec_intersect(double[] pl_point,double[] pl_vec, double[] point, double[] vec )
        {
            bool intersection = false;

            double[] b = new double[2];
            for (int q2 = 0; q2 < 2; q2++)
            { b[q2] = point[q2] - pl_point[q2]; }

            double[,] A = new double[2, 2] { { -vec[0],pl_vec[0] }, { -vec[1], pl_vec[1] } };

            double[] lamda = invert2by2(A, b);

            if (lamda[0] <= 1)
            {
                if (lamda[0] >= 0)
                {
                    if (lamda[1] <= 1)
                    {
                        if (lamda[1] >= 0)
                        { intersection = true; }
                    }
                }
            }

            return intersection;

        }

        public static double [] invert2by2 (double [,] A, double [] b )
        {
            double[,] A_inverse = new double[2, 2] { {A[1,1],-A[0,1] }, {-A[1,0], A[0,0] } };
            for (int q2 = 0; q2 < 2; q2++)
            {
                for (int q1 = 0; q1 < 2; q1++)
                { A_inverse[q1, q2] = A_inverse[q1, q2] / (A[0,0]*A[1,1]-A[0,1]*A[1,0]); }
            }

            double[] x = RVEExamplesBuilder.MatVecMult(A_inverse, b);
            return x;

        }

        public static double[] invert3by3(double[,] A, double[] b)
        {
            double[,] A_inverse = new double[3, 3] { { A[2, 2]*A[1,1]- A[2, 1] * A[1, 2],-(A[2, 2] * A[0, 1] - A[2, 1] * A[0, 2]), A[1,2] * A[0, 1] - A[1, 1] * A[0, 2] },
                { -(A[2,2]*A[1,0]-A[2,0]*A[1,2]),A[2,2]*A[0,0]-A[2,0]*A[0,2], -(A[1,2]*A[0,0]-A[1,0]*A[0,2]) },
                {A[2,1]*A[1,0]-A[2,0]*A[1,1],-(A[2,1]*A[0,0]-A[2,0]*A[0,1]),A[1,1]*A[0,0]-A[1,0]*A[0,1] } };



            double det_A;

            det_A = A[0, 0] * (A[2, 2] * A[1, 1] - A[2, 1] * A[1, 2]) - A[1, 0] * (A[2, 2] * A[0, 1] - A[2, 1] * A[0, 2]) + A[2, 0] * (A[1, 2] * A[0, 1] - A[1, 1] * A[0, 2]);
            for (int q2 = 0; q2 < 3; q2++)
            {
                for (int q1 = 0; q1 < 3; q1++)
                { A_inverse[q1, q2] = A_inverse[q1, q2] / det_A; }
            }

            double[] x = RVEExamplesBuilder.MatVecMult(A_inverse, b);
            return x;

        }

        public static bool combine_bool_any(bool [,] matrix)
        {
            bool output_value = false;
            bool stop_outer_loop = false;
            for (int i1 = 0; i1 < matrix.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < matrix.GetLength(1); i2++)
                {
                    output_value = (output_value || matrix[i1, i2]);
                    if (output_value)
                    {
                        stop_outer_loop = true;
                        break;
                    }
                }
                if (stop_outer_loop)
                {
                    break;
                }
            }

            return output_value;
        }

        public static bool combine_bool_any(bool[] array)
        {
            bool output_value = false;           
            for (int i1 = 0; i1 < array.GetLength(0); i1++)
            {
                output_value = (output_value || array[i1]);
                if (output_value)
                {
                    break;
                }
            }

            return output_value;
        }

        public static double sum_of_bool_array_elements(bool [,] array)
        {
            double output_value = 0;
            
            for (int i1 = 0; i1 < array.GetLength(0); i1++)
            {
                for (int i2 =0; i2 < array.GetLength(1); i2++)
                {
                    if (array[i1, i2])
                    { output_value += 1; }
                    
                }
            }

            return output_value;
        }

        public static double[,] Clone_array(double [,] matrix)
        {
            double[,] array_copy = new double[matrix.GetLength(0), matrix.GetLength(1)];
            for (int i1 = 0; i1 < matrix.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < matrix.GetLength(1); i2++)
                {
                    array_copy[i1, i2] = matrix[i1, i2];
                }
            }
            return array_copy;
        }

        public static bool[,] Clone_array(bool[,] matrix)
        {
            bool[,] array_copy = new bool[matrix.GetLength(0), matrix.GetLength(1)];
            for (int i1 = 0; i1 < matrix.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < matrix.GetLength(1); i2++)
                {
                    array_copy[i1, i2] = matrix[i1, i2];
                }
            }
            return array_copy;
        }
        //public static 

        public static double rand()
        {
            double random1=0.5;
            return random1;
        }

        //copy edw dhmiourgia tous 25_000 kai tous 50_000
        public static void Reference2RVEExample50_000withRenumberingwithInput(Model model)
        {
            double[,] Dq;
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_50_000_renu_new_multiple_algorithms_check_develop\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_50_000_renu_new_multiple_algorithms_check_develop\Fxk_p_komvoi_rve.txt";
            string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_50_000_renu_new_multiple_algorithms_check_develop\o_xsunol_gs_{0}.txt";
            int subdiscr1 = 6;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 23;
            int subdiscr1_shell = 14;
            int discr1_shell = 1;
            mpgp = RVEkanoninkhsGewmetriasBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;


            int graphene_sheets_number = 10;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
            double[][] ekk_xyz = new double[graphene_sheets_number][];
            //double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            RVEkanoninkhsGewmetriasBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);

            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample100_000withRenumberingwithInput(Model model)
        {
            double[,] Dq;
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_100_000_renu_new_multiple_algorithms_check_develop\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_100_000_renu_new_multiple_algorithms_check_develop\Fxk_p_komvoi_rve.txt";
            string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_100_000_renu_new_multiple_algorithms_check_develop\o_xsunol_gs_{0}.txt";
            int subdiscr1 = 8;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 28;
            int subdiscr1_shell = 18;
            int discr1_shell = 1;
            mpgp = RVEkanoninkhsGewmetriasBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;


            int graphene_sheets_number = 10;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
            double[][] ekk_xyz = new double[graphene_sheets_number][];
            //double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            RVEkanoninkhsGewmetriasBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);

            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

    }
}
