﻿using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Elements.Boundary;
using ISAAR.MSolve.IGA.Elements.Structural;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Geometry;
using ISAAR.MSolve.IGA.Geometry.NurbsMesh;
using ISAAR.MSolve.IGA.Loading.LineLoads;
using ISAAR.MSolve.IGA.Loading.LoadElementFactories;
using ISAAR.MSolve.IGA.Loading.NodalLoads;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interpolation;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using Newtonsoft.Json;
using Xunit;
using MatlabWriter = ISAAR.MSolve.LinearAlgebra.Output.MatlabWriter;

namespace ISAAR.MSolve.IGA.Tests
{
    public class NurbsKirchhoffLoveShells
	{
		private List<ControlPoint> ElementControlPoints()
		{
			return new List<ControlPoint>
			{
				new ControlPoint {ID = 0, X = 0.0, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 1, X = 0.0, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 2, X = 0.0, Y =  1.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 3, X = 16.66666667, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 4, X = 16.66666667, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 5, X = 16.66666667, Y =  1.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 6, X = 33.33333333, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 7, X = 33.33333333, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 8, X = 33.33333333, Y =  1.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 9, X = 50.0, Y =  0.0, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 10, X = 50.0, Y =  0.5, Z = 0.0, WeightFactor = 1.0},
				new ControlPoint {ID = 11, X = 50.0, Y =  1.0, Z = 0.0, WeightFactor = 1.0},
			};
		}

		private List<Knot> ElementKnot()
		{
			return new List<Knot>()
			{
				new Knot() {ID = 0, Ksi = 0.0, Heta = 0.0, Zeta = 0.0},
				new Knot() {ID = 1, Ksi = 0.0, Heta = 1.0, Zeta = 0.0},
				new Knot() {ID = 2, Ksi = 1.0, Heta = 0.0, Zeta = 0.0},
				new Knot() {ID = 3, Ksi = 1.0, Heta = 1.0, Zeta = 0.0}
			};
		}

		private double[] KnotValueVectorKsi()
		{
			return new double[8]
			{
				0, 0, 0, 0, 1, 1, 1, 1
			};
		}

		private double[] KnotValueVectorHeta()
		{
			return new double[6]
			{
				0, 0, 0, 1, 1, 1
			};
		}

		private KirchhoffLoveShell Element
		{
			get
			{
                var thickness = 1;
                var degreeKsi = 3;
                var degreeHeta = 2;
                var numberOfControlPointsHeta = 3;
                var knotValueVectorKsi = KnotValueVectorKsi();
                var knotValueVectorHeta = KnotValueVectorHeta();

				var gauss = new GaussQuadrature();
                var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, ElementKnot()).ToArray();
				var nurbs = new Nurbs2D(degreeKsi, knotValueVectorKsi, degreeHeta, knotValueVectorHeta,
                    ElementControlPoints().ToArray(), gaussPoints);
                var material = new ShellElasticSectionMaterial2D()
                {
                    YoungModulus = 100,
                    PoissonRatio = 0.0
                };
                var element = new KirchhoffLoveShell(material, nurbs, gaussPoints, thickness);
				var patch = new Patch();
				
				foreach (var controlPoint in ElementControlPoints())
					element.ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
				foreach (var knot in ElementKnot())
					element.KnotsDictionary.Add(knot.ID, knot);
				
				element.Patch = patch;
				return element;
			}
		}

		private readonly Matrix _expectedStiffnessMatrix = Matrix.CreateFromArray(new double[36, 36]
		{
			{
				476.9104762300273, 12.500000000000012, 0.0, -237.7352381150549, 8.333333333333337, 0.0,
				-237.97523811504095, 4.166666666666667, 0.0, 237.73523809556187, -7.500000000000006, 0.0,
				-119.2276190477295, -4.999999999999999, 0.0, -119.10761904774655, -2.5000000000000004, 0.0,
				94.99809523305787, -3.750000000000005, 0.0, -47.73904761653239, -2.5000000000000013, 0.0,
				-47.6590476165312, -1.2500000000000007, 0.0, 23.68952380802038, -1.250000000000001, 0.0,
				-11.964761904017053, -0.8333333333333335, 0.0, -11.924761904014751, -0.4166666666666668, 0.0
			},
			{
				12.50000000000001, 952.740952460116, 0.0, -8.333333333333346, -476.010476230079, 0.0,
				-4.166666666666668, -476.1304762300716, 0.0, 7.500000000000005, 476.0104761910465, 0.0,
				-5.000000000000006, -238.18523809549762, 0.0, -2.500000000000001, -238.12523809550598, 0.0,
				3.7500000000000036, 190.35619046612086, 0.0, -2.5000000000000036, -95.2980952330622, 0.0,
				-1.2500000000000004, -95.25809523306152, 0.0, 1.250000000000001, 47.559047616051025, 0.0,
				-0.8333333333333345, -23.83952380802895, 0.0, -0.41666666666666663, -23.819523808027782, 0.0
			},
			{
				0.0, 0.0, 238.8953981149884, 0.0, 0.0, -476.5903962300457, 0.0, 0.0, 237.69526478172372, 0.0, 0.0,
				118.64737904781246, 0.0, 0.0, -237.89535809553882, 0.0, 0.0, 119.24757904772663, 0.0, 0.0,
				47.35238094986011, 0.0, 0.0, -95.10476189972607, 0.0, 0.0, 47.75238094986587, 0.0, 0.0,
				11.77150857067266, 0.0, 0.0, -23.74281714135675, 0.0, 0.0, 11.971441904017421
			},
			{
				-237.7352381150549, -8.333333333333346, 0.0, 476.67047623004134, -3.400108834891377E-16, 0.0,
				-237.7352381150549, 8.333333333333346, 0.0, -119.2276190477295, 5.000000000000005, 0.0,
				237.85523809554485, 1.3170345890267665E-16, 0.0, -119.22761904772952, -5.000000000000006, 0.0,
				-47.739047616532375, 2.500000000000003, 0.0, 95.07809523305905, 1.3964523981613297E-16, 0.0,
				-47.739047616532396, -2.500000000000003, 0.0, -11.96476190401705, 0.8333333333333344, 0.0,
				23.72952380802267, 1.3986208025063007E-17, 0.0, -11.96476190401705, -0.8333333333333344, 0.0
			},
			{
				8.333333333333337, -476.01047623007895, 0.0, 6.626592856254568E-16, 952.6209524601236, 0.0,
				-8.333333333333336, -476.01047623007895, 0.0, 5.000000000000002, -238.18523809549762, 0.0,
				4.023474262093707E-16, 476.07047619103815, 0.0, -5.000000000000003, -238.18523809549768, 0.0,
				2.5000000000000013, -95.29809523306223, 0.0, 8.304988641238964E-17, 190.39619046612162, 0.0,
				-2.500000000000001, -95.29809523306223, 0.0, 0.8333333333333333, -23.839523808028954, 0.0,
				3.469446951953614E-17, 47.5790476160522, 0.0, -0.8333333333333334, -23.839523808028957, 0.0
			},
			{
				0.0, 0.0, -476.5903962300457, 0.0, 0.0, 953.1810591267578, 0.0, 0.0, -476.5903962300457, 0.0, 0.0,
				-237.89535809553885, 0.0, 0.0, 475.790316191078, 0.0, 0.0, -237.89535809553885, 0.0, 0.0,
				-95.10476189972609, 0.0, 0.0, 190.209523799452, 0.0, 0.0, -95.10476189972609, 0.0, 0.0,
				-23.74281714135675, 0.0, 0.0, 47.48576761604683, 0.0, 0.0, -23.74281714135675
			},
			{
				-237.97523811504095, -4.166666666666669, 0.0, -237.7352381150549, -8.333333333333336, 0.0,
				476.91047623002726, -12.500000000000009, 0.0, -119.10761904774657, 2.4999999999999996, 0.0,
				-119.22761904772953, 4.999999999999998, 0.0, 237.73523809556184, 7.500000000000007, 0.0,
				-47.6590476165312, 1.2500000000000004, 0.0, -47.7390476165324, 2.5000000000000004, 0.0,
				94.99809523305788, 3.7500000000000044, 0.0, -11.924761904014753, 0.41666666666666674, 0.0,
				-11.964761904017053, 0.8333333333333337, 0.0, 23.68952380802038, 1.250000000000001, 0.0
			},
			{
				4.166666666666667, -476.1304762300716, 0.0, 8.333333333333348, -476.01047623007895, 0.0,
				-12.500000000000009, 952.7409524601159, 0.0, 2.5000000000000004, -238.12523809550598, 0.0,
				5.000000000000006, -238.18523809549762, 0.0, -7.500000000000005, 476.0104761910465, 0.0,
				1.2500000000000004, -95.25809523306152, 0.0, 2.500000000000003, -95.2980952330622, 0.0,
				-3.750000000000003, 190.35619046612086, 0.0, 0.4166666666666666, -23.819523808027782, 0.0,
				0.8333333333333345, -23.839523808028957, 0.0, -1.250000000000001, 47.559047616051025, 0.0
			},
			{
				0.0, 0.0, 237.69526478172372, 0.0, 0.0, -476.5903962300457, 0.0, 0.0, 238.89539811498847, 0.0, 0.0,
				119.24757904772665, 0.0, 0.0, -237.89535809553885, 0.0, 0.0, 118.64737904781246, 0.0, 0.0,
				47.75238094986587, 0.0, 0.0, -95.10476189972609, 0.0, 0.0, 47.35238094986011, 0.0, 0.0,
				11.971441904017421, 0.0, 0.0, -23.74281714135675, 0.0, 0.0, 11.771508570672662
			},
			{
				237.73523809556187, 7.500000000000005, 0.0, -119.2276190477295, 5.000000000000002, 0.0,
				-119.10761904774655, 2.5000000000000004, 0.0, 286.19428569914965, 1.721713049906981E-16, 0.0,
				-142.6171428496091, -1.7271340607694086E-16, 0.0, -142.77714284959762, -2.441623292437356E-16, 0.0,
				214.4057142722313, -3.7500000000000036, 0.0, -107.08285713612945, -2.500000000000001, 0.0,
				-107.12285713612478, -1.2500000000000004, 0.0, 94.99809523305787, -3.7500000000000036, 0.0,
				-47.73904761653239, -2.5000000000000004, 0.0, -47.65904761653119, -1.2500000000000002, 0.0
			},
			{
				-7.500000000000006, 476.0104761910465, 0.0, 5.000000000000007, -238.18523809549762, 0.0,
				2.4999999999999996, -238.12523809550598, 0.0, 2.8319360745321376E-16, 571.6685713983505, 0.0,
				-2.96827449772219E-16, -285.59428569919254, 0.0, -2.88560408207017E-16, -285.6742856991866, 0.0,
				3.750000000000003, 428.6314285444833, 0.0, -2.5000000000000036, -214.2557142722486, 0.0, -1.25,
				-214.2757142722461, 0.0, 3.7500000000000036, 190.3561904661209, 0.0, -2.500000000000004,
				-95.29809523306223, 0.0, -1.25, -95.25809523306152, 0.0
			},
			{
				0.0, 0.0, 118.64737904781248, 0.0, 0.0, -237.89535809553888, 0.0, 0.0, 119.24757904772665, 0.0, 0.0,
				143.390956182887, 0.0, 0.0, -285.98071236583155, 0.0, 0.0, 142.59055618294428, 0.0, 0.0,
				107.2759504694407, 0.0, 0.0, -214.35250093890406, 0.0, 0.0, 107.07615046946349, 0.0, 0.0,
				47.35238094986013, 0.0, 0.0, -95.10476189972606, 0.0, 0.0, 47.75238094986587
			},
			{
				-119.22761904772949, -5.000000000000005, 0.0, 237.85523809554485, 4.301030018249996E-16, 0.0,
				-119.22761904772956, 5.000000000000006, 0.0, -142.61714284960914, -7.409708697309059E-16, 0.0,
				286.03428569916116, 2.732189474663471E-16, 0.0, -142.61714284960917, -3.122502256758253E-17, 0.0,
				-107.08285713612948, 2.5000000000000027, 0.0, 214.36571427223606, -3.0184188481996443E-16, 0.0,
				-107.08285713612949, -2.500000000000003, 0.0, -47.73904761653241, 2.5000000000000036, 0.0,
				95.07809523305909, -2.393918396847994E-16, 0.0, -47.73904761653241, -2.5000000000000036, 0.0
			},
			{
				-4.999999999999999, -238.1852380954976, 0.0, 3.537209587733958E-16, 476.07047619103815, 0.0,
				4.999999999999998, -238.18523809549765, 0.0, 4.933119884809045E-17, -285.59428569919254, 0.0,
				3.8424124992886277E-16, 571.5885713983567, 0.0, -4.744468706796567E-16, -285.59428569919265, 0.0,
				2.4999999999999996, -214.2557142722486, 0.0, -1.0408340855860843E-16, 428.6114285444858, 0.0,
				-2.5000000000000004, -214.25571427224864, 0.0, 2.5, -95.29809523306223, 0.0, -2.498001805406602E-16,
				190.39619046612162, 0.0, -2.5, -95.2980952330622, 0.0
			},
			{
				0.0, 0.0, -237.89535809553882, 0.0, 0.0, 475.79031619107803, 0.0, 0.0, -237.89535809553888, 0.0, 0.0,
				-285.98071236583155, 0.0, 0.0, 571.9622247316629, 0.0, 0.0, -285.9807123658316, 0.0, 0.0,
				-214.35250093890406, 0.0, 0.0, 428.70460187780833, 0.0, 0.0, -214.35250093890403, 0.0, 0.0,
				-95.10476189972607, 0.0, 0.0, 190.20952379945214, 0.0, 0.0, -95.10476189972607
			},
			{
				-119.10761904774655, -2.500000000000001, 0.0, -119.22761904772953, -5.000000000000003, 0.0,
				237.73523809556184, -7.500000000000005, 0.0, -142.77714284959762, -2.88560408207017E-16, 0.0,
				-142.61714284960917, -2.5240226575462543E-16, 0.0, 286.19428569914965, 6.453171330633722E-16, 0.0,
				-107.12285713612478, 1.2500000000000002, 0.0, -107.08285713612948, 2.5000000000000013, 0.0,
				214.40571427223134, 3.7500000000000036, 0.0, -47.65904761653119, 1.2500000000000002, 0.0,
				-47.73904761653242, 2.5000000000000004, 0.0, 94.99809523305791, 3.7500000000000036, 0.0
			},
			{
				-2.4999999999999996, -238.12523809550598, 0.0, -5.000000000000006, -238.18523809549768, 0.0,
				7.500000000000005, 476.0104761910465, 0.0, -1.3227266504323154E-16, -285.6742856991866, 0.0,
				-3.122502256758253E-17, -285.59428569919265, 0.0, 4.2327252813834093E-16, 571.6685713983507, 0.0, 1.25,
				-214.2757142722461, 0.0, 2.500000000000004, -214.2557142722486, 0.0, -3.750000000000004,
				428.6314285444833, 0.0, 1.2500000000000002, -95.25809523306152, 0.0, 2.500000000000005,
				-95.29809523306226, 0.0, -3.750000000000004, 190.35619046612092, 0.0
			},
			{
				0.0, 0.0, 119.24757904772665, 0.0, 0.0, -237.89535809553888, 0.0, 0.0, 118.64737904781246, 0.0, 0.0,
				142.59055618294428, 0.0, 0.0, -285.9807123658316, 0.0, 0.0, 143.390956182887, 0.0, 0.0,
				107.07615046946349, 0.0, 0.0, -214.35250093890403, 0.0, 0.0, 107.27595046944072, 0.0, 0.0,
				47.75238094986588, 0.0, 0.0, -95.10476189972611, 0.0, 0.0, 47.352380949860134
			},

			{
				94.99809523305785, 3.7500000000000036, 0.0, -47.739047616532375, 2.500000000000001, 0.0,
				-47.65904761653119, 1.2500000000000004, 0.0, 214.40571427223136, 3.750000000000003, 0.0,
				-107.08285713612945, 2.4999999999999996, 0.0, -107.12285713612475, 1.2499999999999998, 0.0,
				286.1942856991496, -1.3162214373974024E-15, 0.0, -142.6171428496091, -8.396061623727746E-16, 0.0,
				-142.7771428495976, -2.498001805406602E-16, 0.0, 237.73523809556178, -7.500000000000007, 0.0,
				-119.2276190477295, -5.000000000000001, 0.0, -119.10761904774655, -2.5000000000000004, 0.0
			},
			{
				-3.750000000000005, 190.3561904661209, 0.0, 2.500000000000003, -95.29809523306223, 0.0,
				1.2500000000000004, -95.25809523306152, 0.0, -3.7500000000000036, 428.6314285444833, 0.0,
				2.5000000000000027, -214.2557142722486, 0.0, 1.2500000000000002, -214.2757142722461, 0.0,
				-1.2051991349348867E-15, 571.6685713983507, 0.0, -3.3306690738754696E-16, -285.5942856991926, 0.0,
				2.740863092043355E-16, -285.67428569918656, 0.0, 7.500000000000003, 476.01047619104645, 0.0,
				-5.000000000000003, -238.18523809549768, 0.0, -2.4999999999999987, -238.12523809550595, 0.0
			},
			{
				0.0, 0.0, 47.35238094986012, 0.0, 0.0, -95.10476189972609, 0.0, 0.0, 47.752380949865874, 0.0, 0.0,
				107.2759504694407, 0.0, 0.0, -214.35250093890406, 0.0, 0.0, 107.07615046946349, 0.0, 0.0,
				143.39095618288698, 0.0, 0.0, -285.98071236583155, 0.0, 0.0, 142.5905561829443, 0.0, 0.0,
				118.64737904781246, 0.0, 0.0, -237.89535809553877, 0.0, 0.0, 119.24757904772659
			},
			{
				-47.73904761653238, -2.500000000000003, 0.0, 95.07809523305906, 8.304988641238964E-17, 0.0,
				-47.739047616532396, 2.500000000000003, 0.0, -107.08285713612945, -2.5000000000000036, 0.0,
				214.36571427223606, -2.220446049250313E-16, 0.0, -107.08285713612948, 2.500000000000004, 0.0,
				-142.6171428496091, -1.1102230246251565E-16, 0.0, 286.0342856991613, -9.43689570931383E-16, 0.0,
				-142.61714284960917, 8.881784197001252E-16, 0.0, -119.22761904772956, 5.000000000000007, 0.0,
				237.85523809554488, -8.326672684688674E-16, 0.0, -119.22761904772955, -5.000000000000007, 0.0
			},
			{
				-2.500000000000001, -95.2980952330622, 0.0, 1.3997050046787862E-16, 190.39619046612162, 0.0,
				2.5000000000000004, -95.29809523306223, 0.0, -2.5000000000000004, -214.2557142722486, 0.0,
				-1.8735013540549517E-16, 428.6114285444858, 0.0, 2.5000000000000004, -214.2557142722486, 0.0,
				-1.061650767297806E-15, -285.5942856991926, 0.0, -9.992007221626409E-16, 571.588571398357, 0.0,
				9.992007221626409E-16, -285.5942856991927, 0.0, 4.999999999999997, -238.18523809549768, 0.0,
				-4.440892098500626E-16, 476.0704761910384, 0.0, -4.999999999999996, -238.1852380954977, 0.0
			},
			{
				0.0, 0.0, -95.10476189972609, 0.0, 0.0, 190.20952379945203, 0.0, 0.0, -95.10476189972609, 0.0, 0.0,
				-214.35250093890406, 0.0, 0.0, 428.70460187780833, 0.0, 0.0, -214.35250093890403, 0.0, 0.0,
				-285.98071236583155, 0.0, 0.0, 571.9622247316629, 0.0, 0.0, -285.9807123658316, 0.0, 0.0,
				-237.89535809553888, 0.0, 0.0, 475.79031619107803, 0.0, 0.0, -237.89535809553885
			},
			{
				-47.65904761653119, -1.2500000000000004, 0.0, -47.739047616532396, -2.500000000000001, 0.0,
				94.99809523305787, -3.750000000000003, 0.0, -107.12285713612476, -1.25, 0.0, -107.08285713612948,
				-2.5000000000000004, 0.0, 214.40571427223134, -3.750000000000004, 0.0, -142.7771428495976,
				2.8102520310824275E-16, 0.0, -142.61714284960917, 7.771561172376096E-16, 0.0, 286.19428569914965,
				-3.552713678800501E-15, 0.0, -119.10761904774655, 2.5, 0.0, -119.22761904772966, 5.000000000000003, 0.0,
				237.7352380955619, 7.50000000000001, 0.0
			},
			{
				-1.2500000000000004, -95.25809523306152, 0.0, -2.5000000000000036, -95.29809523306223, 0.0,
				3.750000000000005, 190.3561904661209, 0.0, -1.2500000000000002, -214.2757142722461, 0.0,
				-2.500000000000003, -214.25571427224864, 0.0, 3.7500000000000044, 428.63142854448336, 0.0,
				-2.498001805406602E-16, -285.67428569918656, 0.0, 8.881784197001252E-16, -285.5942856991927, 0.0,
				-3.1086244689504383E-15, 571.6685713983508, 0.0, 2.499999999999999, -238.12523809550598, 0.0,
				5.000000000000008, -238.18523809549794, 0.0, -7.500000000000003, 476.0104761910467, 0.0
			},
			{
				0.0, 0.0, 47.752380949865874, 0.0, 0.0, -95.10476189972609, 0.0, 0.0, 47.35238094986012, 0.0, 0.0,
				107.07615046946349, 0.0, 0.0, -214.35250093890403, 0.0, 0.0, 107.27595046944072, 0.0, 0.0,
				142.5905561829443, 0.0, 0.0, -285.9807123658316, 0.0, 0.0, 143.39095618288707, 0.0, 0.0,
				119.24757904772672, 0.0, 0.0, -237.895358095539, 0.0, 0.0, 118.64737904781254
			},
			{
				23.689523808020375, 1.250000000000001, 0.0, -11.964761904017053, 0.8333333333333333, 0.0,
				-11.924761904014751, 0.41666666666666663, 0.0, 94.99809523305787, 3.7500000000000036, 0.0,
				-47.739047616532396, 2.5, 0.0, -47.65904761653119, 1.25, 0.0, 237.73523809556178, 7.500000000000003,
				0.0, -119.22761904772956, 4.999999999999997, 0.0, -119.10761904774655, 2.499999999999999, 0.0,
				476.9104762300271, -12.500000000000007, 0.0, -237.73523811505478, -8.33333333333333, 0.0,
				-237.97523811504078, -4.166666666666666, 0.0
			},
			{
				-1.250000000000001, 47.55904761605103, 0.0, 0.8333333333333343, -23.839523808028954, 0.0,
				0.41666666666666674, -23.819523808027782, 0.0, -3.7500000000000036, 190.3561904661209, 0.0,
				2.5000000000000036, -95.29809523306223, 0.0, 1.2500000000000002, -95.25809523306152, 0.0,
				-7.500000000000007, 476.01047619104645, 0.0, 5.000000000000007, -238.18523809549774, 0.0, 2.5,
				-238.12523809550595, 0.0, -12.500000000000007, 952.7409524601159, 0.0, 8.333333333333343,
				-476.0104762300789, 0.0, 4.166666666666665, -476.13047623007134, 0.0
			},
			{
				0.0, 0.0, 11.771508570672664, 0.0, 0.0, -23.742817141356753, 0.0, 0.0, 11.971441904017421, 0.0, 0.0,
				47.35238094986012, 0.0, 0.0, -95.10476189972607, 0.0, 0.0, 47.75238094986588, 0.0, 0.0,
				118.64737904781245, 0.0, 0.0, -237.89535809553888, 0.0, 0.0, 119.24757904772672, 0.0, 0.0,
				238.89539811498832, 0.0, 0.0, -476.59039623004537, 0.0, 0.0, 237.69526478172352
			},
			{
				-11.96476190401705, -0.8333333333333345, 0.0, 23.72952380802267, 6.418476861114186E-17, 0.0,
				-11.964761904017053, 0.8333333333333345, 0.0, -47.73904761653238, -2.500000000000004, 0.0,
				95.07809523305907, -5.551115123125783E-17, 0.0, -47.73904761653241, 2.5000000000000044, 0.0,
				-119.2276190477295, -5.000000000000003, 0.0, 237.85523809554488, -2.220446049250313E-16, 0.0,
				-119.22761904772968, 5.000000000000008, 0.0, -237.73523811505478, 8.333333333333343, 0.0,
				476.6704762300413, -3.552713678800501E-15, 0.0, -237.73523811505493, -8.333333333333343, 0.0
			},
			{
				-0.8333333333333335, -23.839523808028954, 0.0, 3.2526065174565133E-19, 47.5790476160522, 0.0,
				0.8333333333333337, -23.839523808028957, 0.0, -2.5000000000000004, -95.29809523306223, 0.0,
				-2.8796409701215E-16, 190.39619046612162, 0.0, 2.500000000000001, -95.29809523306226, 0.0, -5.0,
				-238.18523809549765, 0.0, -9.992007221626409E-16, 476.0704761910384, 0.0, 5.000000000000003,
				-238.18523809549794, 0.0, -8.33333333333333, -476.01047623007884, 0.0, -3.3306690738754696E-15,
				952.6209524601236, 0.0, 8.333333333333332, -476.01047623007906, 0.0
			},
			{
				0.0, 0.0, -23.742817141356753, 0.0, 0.0, 47.48576761604683, 0.0, 0.0, -23.742817141356753, 0.0, 0.0,
				-95.10476189972606, 0.0, 0.0, 190.20952379945214, 0.0, 0.0, -95.1047618997261, 0.0, 0.0,
				-237.89535809553877, 0.0, 0.0, 475.79031619107803, 0.0, 0.0, -237.895358095539, 0.0, 0.0,
				-476.59039623004537, 0.0, 0.0, 953.181059126757, 0.0, 0.0, -476.5903962300451
			},
			{
				-11.924761904014753, -0.41666666666666663, 0.0, -11.964761904017053, -0.8333333333333334, 0.0,
				23.689523808020375, -1.250000000000001, 0.0, -47.65904761653119, -1.25, 0.0, -47.73904761653239,
				-2.4999999999999996, 0.0, 94.99809523305791, -3.750000000000004, 0.0, -119.10761904774655,
				-2.4999999999999987, 0.0, -119.22761904772955, -4.9999999999999964, 0.0, 237.7352380955619,
				-7.500000000000003, 0.0, -237.97523811504078, 4.166666666666665, 0.0, -237.73523811505493,
				8.333333333333332, 0.0, 476.91047623002703, 12.500000000000005, 0.0
			},
			{
				-0.4166666666666668, -23.819523808027782, 0.0, -0.8333333333333344, -23.839523808028957, 0.0,
				1.250000000000001, 47.559047616051025, 0.0, -1.2500000000000002, -95.25809523306152, 0.0,
				-2.500000000000004, -95.2980952330622, 0.0, 3.7500000000000036, 190.35619046612092, 0.0,
				-2.5000000000000004, -238.12523809550592, 0.0, -5.000000000000007, -238.1852380954977, 0.0,
				7.50000000000001, 476.0104761910467, 0.0, -4.166666666666666, -476.13047623007134, 0.0,
				-8.333333333333343, -476.010476230079, 0.0, 12.500000000000005, 952.7409524601155, 0.0
			},
			{
				0.0, 0.0, 11.971441904017421, 0.0, 0.0, -23.742817141356753, 0.0, 0.0, 11.771508570672662, 0.0, 0.0,
				47.75238094986587, 0.0, 0.0, -95.10476189972609, 0.0, 0.0, 47.352380949860134, 0.0, 0.0,
				119.24757904772656, 0.0, 0.0, -237.89535809553885, 0.0, 0.0, 118.64737904781254, 0.0, 0.0,
				237.69526478172352, 0.0, 0.0, -476.5903962300451, 0.0, 0.0, 238.89539811498815
			}
		});

		private const double Tolerance = 1e-11;

		[Fact]
		private void TestNurbsKLShellsStiffnessMatrix()
		{
			var element = Element;

			var stiffnessMatrix = element.StiffnessMatrix(element);


			for (var i = 0; i < 36; i++)
			{
				for (var j = 0; j < 36; j++)
				{
					Assert.True(Utilities.AreValuesEqual(_expectedStiffnessMatrix[i, j], stiffnessMatrix[i, j],
						Tolerance));
				}
			}
		}

		[Fact]
		public void IsogeometricCantileverShell()
		{
			string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "CantileverShell.txt");
            var material = new ShellElasticSectionMaterial2D()
            {
                YoungModulus = 100,
                PoissonRatio = 0
            };
			var modelReader = new IsogeometricShellReader(GeometricalFormulation.Linear,filename,sectionMaterial: material);
			var model=modelReader.GenerateModelFromFile();

			model.Loads.Add(new NodalLoad(model.ControlPoints.ToList()[9],StructuralDof.TranslationZ,-1));
            model.Loads.Add(new NodalLoad(model.ControlPoints.ToList()[10],StructuralDof.TranslationZ,-1));
            model.Loads.Add(new NodalLoad(model.ControlPoints.ToList()[11],StructuralDof.TranslationZ,-1));

			for (int i = 0; i < 6; i++)
			{
				model.ControlPointsDictionary[i].Constraints.Add(new Constraint() {DOF = StructuralDof.TranslationX});
				model.ControlPointsDictionary[i].Constraints.Add(new Constraint() {DOF = StructuralDof.TranslationY});
				model.ControlPointsDictionary[i].Constraints.Add(new Constraint() {DOF = StructuralDof.TranslationZ});
			}

			// Solvers
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var expectedSolution = new double[18]
			{
				0.0, 0.0, -7499.999986865148, 0.0, 0.0, -7499.99998660616, 0.0, 0.0, -7499.999986347174, 0.0, 0.0,
				-14999.999980230163, 0.0, 0.0, -14999.999980050825, 0.0, 0.0, -14999.999979871487
			};
			for (int i = 0; i < expectedSolution.Length; i++)
				Utilities.AreValuesEqual(expectedSolution[i], solver.LinearSystems[0].Solution[i], 7);
		}

        [Fact]
        public static void ScordelisLoShell()
        {
            var filename = "ScordelisLoShell";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(),"InputFiles", $"{filename}.txt")
                .ToString(CultureInfo.InvariantCulture);
			var material = new ShellElasticMaterial2Dtransformationb()
            {
                YoungModulus = 4.3210e8,
                PoissonRatio = 0.0
            };
            var modelReader = new IsogeometricShellReader(GeometricalFormulation.NonLinear, filepath, material);
            var model = modelReader.GenerateModelFromFile();

			//model.SurfaceLoads.Add(new SurfaceDistributedLoad(-90, StructuralDof.TranslationY));

            // Rigid diaphragm for AB
            for (var i = 0; i < 19; i++)
            {
                model.ControlPointsDictionary[i * 19].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[i * 19].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
            }

            // Symmetry for CD
            for (var i = 0; i < 19; i++)
            {
                model.ControlPointsDictionary[i * 19 + 18].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });

                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[i * 19 + 18], StructuralDof.TranslationX),
                    new NodalDof(model.ControlPointsDictionary[i * 19 + 17], StructuralDof.TranslationX)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[i * 19 + 18], StructuralDof.TranslationY),
                    new NodalDof(model.ControlPointsDictionary[i * 19 + 17], StructuralDof.TranslationY)));
            }

            // Symmetry for AD
            for (var j = 0; j < 19; j++)
            {
                model.ControlPointsDictionary[j].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[j], StructuralDof.TranslationY),
                    new NodalDof(model.ControlPointsDictionary[j + 19], StructuralDof.TranslationY)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[j], StructuralDof.TranslationZ),
                    new NodalDof(model.ControlPointsDictionary[j + 19], StructuralDof.TranslationZ)));
            }

            // Solvers
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            var cp = model.ControlPointsDictionary.Values.Last();
            var dofA = model.GlobalDofOrdering.GlobalFreeDofs[cp, StructuralDof.TranslationY];

            var solution = solver.LinearSystems[0].Solution[dofA];
        }

        [Fact]
		public void IsogeometricSquareShell()
		{
            var material = new ShellElasticMaterial2Dtransformationb()
            {
                YoungModulus = 10000000,
                PoissonRatio = 0.0
            };
            var filename = "SquareShell";
            string filepath = Path.Combine(Directory.GetCurrentDirectory(),"InputFiles",$"{filename}.json");
            var jsonReader = new JsonModelReader(filepath, shellMaterial:material);

            var (geometry, model)=jsonReader.ReadGeometryAndCreateModel();
            var rightEdgeLoads=geometry.NurbsSurfacePatches[0]
                .CreateLoadForEdge(model,NurbsSurfaceEdges.Right, new DistributedLineLoad(0, -100, 0));

            geometry.NurbsSurfacePatches[0].ConstraintDofsOfEdge(model, NurbsSurfaceEdges.Left, new List<IDofType>()
            {
                StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
            });
            geometry.NurbsSurfacePatches[0].ConstraintDofsOfEdge(model, NurbsSurfaceEdges.Right, new List<IDofType>()
            {
                StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
            });
            geometry.NurbsSurfacePatches[0].ConstraintDofsOfEdge(model, NurbsSurfaceEdges.Bottom, new List<IDofType>()
            {
                StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
            });
            geometry.NurbsSurfacePatches[0].ConstraintDofsOfEdge(model, NurbsSurfaceEdges.Top, new List<IDofType>()
            {
                StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
            });

            Matrix<double> loadVector = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "SquareShell.mat"), "LoadVector");


			for (int i = 0; i < loadVector.ColumnCount; i += 3)
			{
				var indexCP = i / 3;
                model.Loads.Add(new NodalLoad(model.ControlPoints.ToList()[indexCP], StructuralDof.TranslationX, loadVector.At(0, i)));
                model.Loads.Add(new NodalLoad(model.ControlPoints.ToList()[indexCP], StructuralDof.TranslationY, loadVector.At(0, i+1)));
                model.Loads.Add(new NodalLoad(model.ControlPoints.ToList()[indexCP], StructuralDof.TranslationZ, loadVector.At(0, i+2)));
			}

            // Solvers
            var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

            //var paraview = new ParaviewNurbsShells(model, solver.LinearSystems[0].Solution, filename);
            //paraview.CreateParaview2DFile();

            Matrix<double> solutionVectorExpected =
				MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "SquareShell.mat"), "SolutionVector");

			for (int i = 0; i < solutionVectorExpected.ColumnCount; i++)
				Assert.True(Utilities.AreValuesEqual(solutionVectorExpected.At(0, i), solver.LinearSystems[0].Solution[i],
					1e-9));
		}


    }
}