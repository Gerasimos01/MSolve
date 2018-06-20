using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.Materials;//using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.FEM.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;//using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
// compa
using ISAAR.MSolve.FEM.Entities;

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

            ElasticMaterial3D_v2 material1 = new ElasticMaterial3D_v2()
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

        //public static void Reference2RVEExample10000withRenumberingwithInput_1GrSh_sparse_v2(Model model)
        //{
        //    double[,] Dq;
        //    Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
        //    rveMatrixParameters mp;
        //    grapheneSheetParameters gp;
        //    string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh\REF_new_total_numbering.txt";
        //    string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh\Fxk_p_komvoi_rve.txt";
        //    string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh\o_xsunol_gs_{0}.txt";
        //    int subdiscr1 = 4;
        //    int discr1 = 4;
        //    // int discr2 dn xrhsimopoieitai
        //    int discr3 = 10;
        //    int subdiscr1_shell = 7;
        //    int discr1_shell = 1;
        //    mpgp = RVEkanoninkhsGewmetriasBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
        //    mp = mpgp.Item1;
        //    gp = mpgp.Item2;


        //    int graphene_sheets_number = 1;
        //    o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
        //    double[][] ekk_xyz = new double[graphene_sheets_number][];
        //    //double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


        //    Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
        //    RVEkanoninkhsGewmetriasBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);

        //    int hexaElementsNumber = model.ElementsDictionary.Count();

        //    IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
        //    List<int> EmbeddedElementsIDs = new List<int>();
        //    int element_counter_after_Adding_sheet;
        //    element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
        //    int shellElementsNumber;

        //    for (int j = 0; j < graphene_sheets_number; j++)
        //    {
        //        string file_no = (j + 1).ToString();
        //        string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
        //        RVEExamplesBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
        //        shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
        //        //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
        //        //embdeddedGroup.Concat(embdeddedGroup_adittion);
        //        for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
        //        {
        //            EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
        //        }
        //        element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

        //    }

        //    // model: add loads
        //    RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
        //    //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
        //    // model: add constraints
        //    RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
        //    //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

        //    int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
        //    IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
        //    var embeddedGrouping = new EmbeddedCohesiveGrouping_v2(model, hostGroup, embdeddedGroup);
        //}

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

        public static void Reference2RVEExample10000withRenumberingwithInput_forMS(Model model, Dictionary<int, Node> boundaryNodes)
        {
            //PROELEFSI public static void Reference2RVEExample10000withRenumberingwithInput(Model model)
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


            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            RVEkanoninkhsGewmetriasBuilder.HexaElementsOnlyRVEwithRenumbering_forMS(model, mp, Dq, renumbering_vector_path, boundaryNodes);

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

        public static void HexaElementsOnlyRVEwithRenumbering_forMS(Model model, rveMatrixParameters mp, double[,] Dq, string renumberingVectorPath, Dictionary<int, Node> boundaryNodes)
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

            ElasticMaterial3D_v2 material1 = new ElasticMaterial3D_v2()
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
                boundaryNodes.Add(komvos, model.NodesDictionary[komvos]);
            }


        }

    }

}

