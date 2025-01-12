(* ::Package:: *)

(* << utilities.wl *)



(* ::Section:: *)
(*Others*)


form[x_] := Rasterize@TraditionalForm@x



(* ::Section:: *)
(*Potential  energy  U*)


uNPlot[U_, kxTemp_, pxTemp_, zTemp_:3, erg_ : 1/2, opts : OptionsPattern[
    {Plot, addLine -> True, addPoint -> False}]] := Module[{},
    epilog = {};
    
    
    If[OptionValue[addLine],
        line1 = Line[{{-zTemp, erg}, {zTemp, erg}}];
        epilog = Join[epilog, {Directive[{Dashed}], line1}]
    ];
    Utemp = U /. {\[Kappa]xN -> Rationalize[kxTemp], pxN -> pxTemp};
    If[OptionValue[addPoint],
        zNs = NSolveValues[# == erg && zN \[Element] Reals, zN]& /@ Utemp;
        point = Point[{#, erg}& /@ #]& /@ zNs;
        epilog = Join[epilog, {PointSize[Medium], point}]
    ];
    Plot[Evaluate[Utemp], {zN, -zTemp, zTemp}, PlotRange -> {Full, {0,
         1}}, Epilog -> epilog, Axes -> False, Frame -> True, PlotLabel -> StringForm[
        "\[Kappa]x: ``, px: ``", N @ kxTemp, N @ pxTemp], FrameLabel -> {"z", "U"}, 
        Evaluate[FilterRules[{opts}, Options[Plot]]]]
]


uPlot[kxTemp_, pxTemp_, zTemp_ : 3] := Module[{},
  line1 = Line[{{-zTemp, 1/2}, {zTemp, 1/2}}];
  lineStyle = {Dashed};
  Plot[U /. {kx -> kxTemp, px -> pxTemp}, {z, -zTemp, zTemp}, 
   PlotRange -> {Full, {0, 1}},
   Epilog -> {Directive[lineStyle], line1},
   Axes -> False, Frame -> True, 
   PlotLabel -> StringForm["kx: `` px: ``", N@kxTemp, N@pxTemp],
   FrameLabel -> {"z", "\!\(\*SubscriptBox[\(U\), \(z\)]\)"}
   ]
  ]


(* ::Section::Closed:: *)
(*Solver*)


getInitialPoints::usage = 
  "Get Initial Points given equations (need to set `x0, px0, z0, pz0`)"

ClearAll[getInitialPoints]
getInitialPoints[opts : OptionsPattern[]] := Block[{freeVariables},
  freeVariables = DeleteDuplicates@Cases[tempEqs, _Symbol, \[Infinity]];
  (*allVariables={x,px,z,pz};*)
  (*fixedVariables=DeleteElements[allVariables,freeVariables];*)
  
  reg = ImplicitRegion[tempEqs, Evaluate@freeVariables] // DiscretizeRegion;
  If[OptionValue[InitialRegionPlot],
   Print@Region[reg, Axes -> True, AxesLabel -> freeVariables]];
  
  RandomPoint[reg, OptionValue[InitialPoints]] //
     Map[AssociationThread[ToString /@ freeVariables, #] &] //
    Map[Join[<|"x" -> x0, "px" -> px0, "z" -> z0, "pz" -> pz0|>, #] &] (* // Values *)
  ]

Options[getInitialPoints] = {
   InitialPoints -> 8,
   InitialRegionPlot -> False
   };


(* ::Subsubsection:: *)
(*getPoint : Solve the integration using initial conditions*)


ClearAll[getPoint]
getPoint::usage = "Solve the integration using initial conditions"

initialBoundaryRatio = 1.000001;

getPoint[initPoint_, opts : OptionsPattern[]] := 
 Module[{sol, points, when, stopWhen, stop = False},
  (*Test for parallelization*)
  (*{px0,z0,pz0}={px,z,pz}/.First@FindInstance[tempEqs,{px,z,pz},Reals];*)
  (*{px0,z0,pz0}=RandomPoint@DiscretizeRegion@reg;*)
  
  x0 = initPoint["x"];
  px0 = initPoint["px"];
  z0 = initPoint["z"];
  pz0 = initPoint["pz"];

  (* pointFormat := {k x[t], z[t], px[t], pz[t], t}; *)
  
  when = {};

  pointFormat := <|"kx" -> k x[t], "z" -> z[t], "px" -> px[t], "pz" -> pz[t], "t"->t|>; 

  zPointWhenEvent = WhenEvent[z[t] == 0, Sow[pointFormat, "zPoints"]];
  xTurningPointWhenEvent = WhenEvent[px[t] == 0, Sow[pointFormat, "xTurningPoints"]];
  uBndPointWhenEvent = WhenEvent[(px[t])^2 + (k x[t])^2 == 1, Sow[pointFormat, "uBndPoints"]]; (* These boundaries are defined by the following condition: one of solutions of the equation U = 1/2 is equal to zero. *)

  If[OptionValue[RecordzPoints], AppendTo[when, zPointWhenEvent]];
  If[OptionValue[RecordxTurningPoints], AppendTo[when, xTurningPointWhenEvent]];
  If[OptionValue[RecorduBndPoints], AppendTo[when, uBndPointWhenEvent]];

  If[OptionValue[StopAtInitialBoundary],
   stopWhen = {
     cc[0] == 0,
     (*WhenEvent[z[t]\[Equal]0,{{cc[t]}->{1},"RemoveEvent"}],*)
     (*WhenEvent[x[t]==initialBoundaryRatio x0&&cc[t]==1,{stop=t,
     "StopIntegration"}]*)
     WhenEvent[x[t] > initialBoundaryRatio x0, {stop = t, Sow[pointFormat, "stopPoints"], "StopIntegration"}]
     };
   when = Join[when, stopWhen]
   ];
  
  (*TODO*)
  If[OptionValue[StopAtFarBoundary],
   stopWhen = WhenEvent[k x[t] > farBoundary, {stop = t, "StopIntegration"}];
   AppendTo[when, stopWhen]
   ];
  
  {sol, points} = Reap[
    NDSolve[
     Join[hamiltonianEquations, initCondition, when],
     {x, z, px, pz},
     {t, tmin, tmax}
     (*If[OptionValue[StopAtInitialBoundary],
     DiscreteVariables->Element[cc,{0,1}]]*)
     ],
    _, Rule
    ];
  
  Join[
   <|"stop" -> stop, "initialPoint" -> initPoint, "stopPoints" -> {}, "zPoints" -> {}, 
    "uBndPoints" -> {}, "xTurningPoints" -> {}|>,
   Association@points,
   If[OptionValue[ReturnSolution],
    <|"solution" -> sol|>,
    <||>
    ]
   ]
  ]

Options[getPoint] = {
   ReturnSolution -> False,
   StopAtInitialBoundary -> True,
   StopAtFarBoundary -> False,
   RecordzPoints -> True,
   RecordxTurningPoints -> True,
   RecorduBndPoints -> True
   };


ClearAll[getPoints]
getPoints::usage = "Get points, given initial conditions (optionally)"

getPoints[initPoints_List, opts : OptionsPattern[]] := 
 Block[{points, pointsDS},
  
  points = Parallelize@Map[getPoint, initPoints];

  pointsDS = Dataset@points;
  pointsDS = pointsDS[All, Append[#, {
    "numZPoints" -> Length[#zPoints],
    "numUBndPoints" -> Length[#uBndPoints],
    "numXTurningPoints" -> Length[#xTurningPoints],
    "uTime" -> uTime[#uBndPoints]
    }] &];
  Echo[pointsDS[All, "numUBndPoints"] // DeleteDuplicates // Normal // Sort, "numUBndPoints"];
  Echo[pointsDS[All, "numXTurningPoints"] // DeleteDuplicates // Normal// Sort, "numXTurningPoints"];

  If[OptionValue[ReturnDataset],
   pointsDS,
   points
   ]
  ]

getPoints[opts : OptionsPattern[]] := Block[{initPoints},
  initPoints = getInitialPoints[];
  getPoints[initPoints]
  ]

Options[getPoints] = {
   ReturnDataset -> True
   };


uTime[uBndPoints_] := 
 Switch[Length@uBndPoints, 0, Infinity, 1, 
  First[uBndPoints][[
   5]], _, (Last[uBndPoints] - First[uBndPoints])[[5]]]


(*Using `RandomPoint` directly without discretization would not work:
The specified region appears to be unbounded.Appropriate bounds will be \
automatically computed.Explicit bounds may be specified as a third argument*)

ClearAll[regArea]
trappedStandard = 2;
regArea[opts : OptionsPattern[{regArea, getPoints}]] := 
 Module[{allPoints, transientPoints, trappedPoints, allReg, transientReg, 
   trappedReg},
  
  allPointsDS = getPoints[ReturnDataset -> True];
  
  allPoints = allPointsDS // Part[#, All, "zPoints"] & // Join @@ # &;
  transientPoints = 
   allPointsDS // Select[#numUBndPoints <= trappedStandard &] // 
     Part[#, All, "zPoints"] & // Join @@ # &;
  trappedPoints = 
   allPointsDS // Select[#numUBndPoints > trappedStandard &] // 
     Part[#, All, "zPoints"] & // Join @@ # &;
  
  {allReg, transientReg, trappedReg} = If[Length@# > 0,
      ConcaveHullMesh[#[[All, {1, 3}]] // Normal],
      Null] & /@ {allPoints, transientPoints, trappedPoints};
  
  Association[
   "\[Alpha]" -> \[Alpha], "k" -> k, "x0" -> x0, "kx" -> k x0,
   "allReg" -> Rasterize@allReg, "areaOfAllReg" -> Area@allReg,
   "transientReg" -> Rasterize@transientReg, 
   "areaOfTransientReg" -> Area@transientReg,
   "trappedReg" -> Rasterize@trappedReg, 
   "areaOfTrappedReg" -> Area@trappedReg
   ]
  ]

nDataset[allPointsDS_] := Module[{ds},
  ds = Table[
      Module[{pointsDS, points, avgUTime, reg},
       pointsDS = allPointsDS // Select[#numUBndPoints == num &];
       avgUTime = pointsDS // Part[#, All, "uTime"] & // Mean;
       
       points = pointsDS // Part[#, All, "zPoints"] & // Join @@ # &;
       
       reg = If[Length@points != 0,
         ConcaveHullMesh[points[[All, {1, 3}]] // Normal]];
       
       <|"numUBndPoints" -> num, "reg" -> Rasterize@reg, 
        "area" -> Area@reg, "avgUTime" -> avgUTime|>
       ],
      
      {num, 
       allPointsDS[All, "numUBndPoints"] // DeleteDuplicates // Normal}
      ] // Dataset // Sort;
  
  ds[All, Append[#, {
      "avgAreaByNum" -> 
       If[#numUBndPoints != 0, #area/#numUBndPoints, Null],
      "avgAreaByUTime" -> If[#avgUTime != 0, #area/#avgUTime, Null],
      "avgUTimeByNum" -> 
       If[#numUBndPoints != 0, #avgUTime/#numUBndPoints, Null]
      }] &]
  ]



(* ::Section:: *)
(*Plotting utilities*)


uCurve[zc0_ : 5] := Module[{},
  tempSol = NSolveValues[
    D[U, z] == 0 && rawH == 1/2 /. {z -> zc, pz -> 0, kx -> kxs, px -> pxs}, 
    {kxs, pxs}, Reals
    ];

  ParametricPlot[tempSol, {zc, -zc0, zc0}, 
    PlotLabel -> StringForm["c1: `` c2: ``", N@c1, N@c2],
    Axes -> False, Frame -> True, FrameLabel -> {"\[Kappa] x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)"}]
]

uCurve[aTemp_, kmTemp_, zc0_ : 5] := Block[{\[Alpha], c1, c2},
  \[Alpha] = aTemp;
  c1 = Sqrt[kmTemp^2 + 1];
  c2 = 1/c1;
  uCurve[zc0]
  ]

uCurveContourPlot[] := Module[{},
  ContourPlot[
    Length@NSolveValues[
      rawH == erg /. {kx -> kxTemp, px -> pxTemp, pz -> 0}, z, Reals
      ], 
    {kxTemp, -1, 6}, {pxTemp, -3, 3},
    Contours -> {0, 1, 2, 3, 4}, 
    ContourShading -> {None, None, None, Gray, None, Yellow}, 
    MaxRecursion -> 3,
    FrameLabel -> {"\[Kappa] x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)"}]
]

uncertaintyCurve[] := Block[{zcCond},
  tempSol = 
   NSolveValues[
    D[U, z] == 0 && D[U, {z, 2}] < 0 && rawH == 1/2 /. {z -> zc, 
      pz -> 0, kx -> kxs, px -> pxs}, {kxs, pxs}, Reals];
  If[Length@tempSol > 0,
   zcCond = First[First@tempSol][[2]];
   ParametricPlot[tempSol, {zc, -Last@zcCond, Last@zcCond}, 
    PlotStyle -> Directive[Dashed, Black]]
   ]
  ]

uncertaintyCurve[aTemp_, kmTemp_] := Block[{\[Alpha], c1, c2},
  \[Alpha] = aTemp;
  c1 = Sqrt[kmTemp^2 + 1];
  c2 = 1/c1;
  uncertaintyCurve[]
  ]

arcLen[aTemp_, kmTemp_] := Block[{\[Alpha], c1, c2},
  \[Alpha] = aTemp;
  c1 = Sqrt[kmTemp^2 + 1];
  c2 = 1/c1;
  tempSol = 
   NSolveValues[
    D[U, z] == 0 && D[U, {z, 2}] < 0 && rawH == 1/2 /. {z -> zc, 
      pz -> 0}, {kx, px}, Reals];
  If[Length@tempSol > 0,
   tempSol = First@tempSol;
   zcCond = First[tempSol][[2]];
   2 ArcLength[tempSol, {zc, First@zcCond, Last@zcCond}]
   ]
  ]


ClearAll[overviewPlot];
overviewPlot[s_, opts : OptionsPattern[]] := Block[{variable, plotList},

  variable = Keys@First@Flatten@s;
  tRange = Flatten[{t, variable["Domain"] /. s}];

  plot1 = 
   Plot[Evaluate[{k x[t], px[t]} /. s], tRange, 
    PlotLegends -> {"\!\(\*SubscriptBox[\(\[Kappa]\), \(n\)]\) x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)"}, AxesLabel -> Automatic, 
    PlotRange -> Full];
  plot2 = 
   Plot[Evaluate[{z[t], pz[t]} /. s], tRange, 
    PlotLegends -> {"z", "pz"}, AxesLabel -> Automatic, 
    PlotRange -> Full];
  (* plot3 = ResourceFunction["DirectionParametricPlot"][
    Evaluate[{k x[t], px[t]} /. s], Evaluate@tRange, 
    AxesLabel -> {"\!\(\*SubscriptBox[\(\[Kappa]\), \(n\)]\) x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)"}, PlotRange -> Full,
    ColorFunction -> Function[{x, y, t}, GrayLevel[t]],
    Epilog -> Circle[{0, 0}],
    "ArrowSize" -> Medium
    ]; *)

  plot3 = ParametricPlot[
    Evaluate[{k x[t], px[t]} /. s], Evaluate@tRange, 
    AxesLabel -> {"\!\(\*SubscriptBox[\(\[Kappa]\), \(n\)]\) x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)"}, PlotRange -> Full,
    ColorFunction -> Function[{x, y, t}, GrayLevel[t]],
    Epilog -> {Circle[{0, 0}], Line[{{k x0, -10}, {k x0, 10}}]}
    ];

  If[OptionValue[PlotUncertaintyCurve],
    plot3 = Show[plot3, uncertaintyCurve[]];
    ];

  
  plotIz := ListPlot[izt, 
    PlotLegends -> {"Iz"}, 
    AxesLabel -> Automatic, PlotRange -> All];

  kxpxPlot := 
  ParametricPlot[Evaluate[{k x[t], px[t]} /. s], Evaluate@tRange,
   AxesLabel -> {"\!\(\*SubscriptBox[\(\[Kappa]\), \(n\)]\) x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)"},
   ColorFunction -> Function[{x, y, t}, ColorData["Rainbow"][izt@t]], ColorFunctionScaling -> False
   ];
  kxpxIzPlot := Legended[Show[kxpxPlot, uncertaintyCurve[]], BarLegend["Rainbow", LegendLabel -> "Iz"]];

  vy := k x[t] - \[Alpha]/2 z[t]^2;
  plotvy := 
   ParametricPlot[Evaluate[{z[t], vy} /. s], Evaluate@tRange, 
    AxesLabel -> {"z", "vy"}];
  (* plotvy := Plot[vy /. sol, tRange, AxesLabel -> {"t", "vy"}]; *)

  plot3D := 
   ParametricPlot3D[Evaluate[{k x[t], px[t], z[t]} /. sol], 
    Evaluate@tRange,
    AxesLabel -> {"\[Kappa] x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)", 
      "z"}];
  
  plotList = {plot1, plot2, plot3};

  If[OptionValue[PlotIz] || OptionValue[PlotkxpxIz],
    izt = izTInterpolation[s];
  ];

  If[OptionValue[PlotIz],
    AppendTo[plotList, plotIz];
  ];

  If[OptionValue[PlotkxpxIz],
    kxpxIzPlot // Print;
    ];

  If[OptionValue[PlotVy],
    AppendTo[plotList, plotvy];
    ];

  If[OptionValue[Plot3D],
    AppendTo[plotList, plot3D];
    ];
  
  If[OptionValue[Graphis],
    GraphicsRow[Rasterize /@ plotList, ImageSize -> Scaled[1/GoldenRatio]],
    plotList
    ]
  (*ListPlot[ifn] and ListLinePlot[
  ifn] will also plot an InterpolatingFunction ifn directly*)
  (*ListPlot[x/.s]*)
  ]

Options[overviewPlot] = {
  PlotUncertaintyCurve -> False,
  PlotIz -> False,
  PlotkxpxIz -> False,
  PlotVy -> False,
  Plot3D -> False,
  Graphis -> True
};


(* izPlot := Module[{tpoints, tRange, tmin, tmax, tt, ergT, Iz, td},
   tpoints = 200;
   	
   tRange := Flatten[{t, x["Domain"] /. s}];
   tmin = tRange[[2]];
   tmax = tRange[[3]];
   tt = Range[tmin, tmax, (tmax - tmin)/tpoints];
   
   
   ergT = 
    First@(H /. {x[t] -> (x[#] /. s), z[t] -> (z[#] /. s), 
          px[t] -> (px[#] /. s), pz[t] -> (pz[#] /. s)}) & /@ tt;
   Iz = Area[
       ImplicitRegion[
        First[H /. {x[t] -> (x[#] /. s), px[t] -> (px[#] /. s)}] <= 
         erg, {z[t], pz[t]}
        ]
       ] & /@ tt;
   
   td = TemporalData[{ergT, Iz}, {tt}];
   ListPlot[td, 
    PlotLegends -> {"Erg", "\!\(\*SubscriptBox[\(I\), \(z\)]\)"}, 
    AxesLabel -> Automatic, PlotRange -> All]
   ]; *)

ClearAll[plotArea]
plotArea::usage= "Plot Regular, Trapped And Transient Area"

plotArea[allPointsDS_, opts : OptionsPattern[]] := 
  Block[{plotList},
 (* Block[{trappedPoints, transientPoints, regularPoints, zPoints, 
   nonStopPoints, stopPoints, zzPoints, plotList}, *)
  
  trappedPoints := 
   allPointsDS // 
      Select[#numUBndPoints > OptionValue[TrappedStandard] &] // 
     Part[#, All, "zPoints"] & // Apply[Join];
  transientPoints := 
   allPointsDS // 
      Select[0 < #numUBndPoints <= 
         OptionValue[TrappedStandard] &] // 
     Part[#, All, "zPoints"] & // Apply[Join];
  regularPoints := 
   allPointsDS // Select[#numUBndPoints == 0 &] // 
     Part[#, All, "zPoints"] & // Apply[Join];
  zPoints = 
   Select[<|"Regular" -> regularPoints[[All, {1, 3}]], 
     "Trapped" -> trappedPoints[[All, {1, 3}]], 
     "Transient" -> transientPoints[[All, {1, 3}]]|>, 
    Length@# > 0 &];
  
  nonStopPoints := 
   allPointsDS // Select[#stop == False &] // 
     Part[#, All, "zPoints"] & // Apply[Join];
  stopPoints := 
   allPointsDS // Select[#stop > 0 &] // Part[#, All, "zPoints"] & // 
    Apply[Join];
  zzPoints = 
   Select[<|"Stopped" -> stopPoints[[All, {1, 3}]], 
     "Non-Stopped" -> nonStopPoints[[All, {1, 3}]]|>, 
    Length@# > 0 &];
  
  plot01 := ListPlot[zPoints,
    PlotLabel -> "Categorized by num of crossing across the U boundary",
    AxesLabel -> {"\[Kappa] x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)"},
     AspectRatio -> 1];

  plot02 := ListPlot[zzPoints,
    PlotLabel -> "Categorized by whether crossing the current sheet boundary",
    AxesLabel -> {"\[Kappa] x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)"},
     AspectRatio -> 1];
  
  plot03 := ListPlot[
    allPointsDS[All, {"numUBndPoints" -> ToString}][
      GroupBy["numUBndPoints"], Apply[Join], "zPoints", 
      Part[#, All, {1, 3}] &] // Select[Length@# > 0 &],
    AxesLabel -> {"\[Kappa] x", "\!\(\*SubscriptBox[\(p\), \(x\)]\)"},
     AspectRatio -> 1,
    PlotLabel -> "Categorized by num of crossing across the U boundary",
    PlotLegends -> 
     PointLegend[Automatic, LegendLabel -> "Num of uncertain points"]
    ];
  
  plotList = {};
  If[Length@zPoints != 0,
   plot01 // AppendTo[plotList, #] &
   ];
  
  If[Length@zzPoints != 0,
   plot02 // AppendTo[plotList, #] &
   ];
  
  If[OptionValue[Plot03],
   plot03 // Rasterize // Print
   ];
  
  If[Length@nonStopPoints == 0,
   Echo["No non-stop particles"]
   ];
  If[Length@stopPoints == 0,
   Echo["No stop particles"]
   ];
  If[Length@transientPoints == 0,
   Echo["No transient particles"]
   ];
  If[Length@trappedPoints == 0,
   Echo["No trapped particles"]
   ];
  
  If[OptionValue[Graphis],
    GraphicsRow[Rasterize /@ plotList, ImageSize -> Scaled[1/GoldenRatio]],
    plotList
    ]
  ]

Options[plotArea] = {
   TrappedStandard -> 2,
   Plot03 -> False,
   Graphis -> False
   };

(* TODO: izT is not reliable when multiple potential well exists *)
izT[t_, s_] := 
  Area@ImplicitRegion[
    First[rawHH /. {x -> (x[t] /. s), px -> (px[t] /. s)}] <= 
      erg, {z, pz}]/ (2 \[Pi]);

izTInterpolation[s_, opts : OptionsPattern[]] := 
 Module[{variable, tRange, tMin, tMax},
  tpoints = OptionValue[tPoints];
  variable = Keys@First@Flatten@s;
  tRange = Flatten[{t, variable["Domain"] /. s}];
  tMin = tRange[[2]];
  tMax = tRange[[3]];
  ts = Range[tMin, tMax, (tMax - tMin)/tpoints];
  izTs = ParallelMap[izT[#, s] &, ts];
  Switch[OptionValue[ReturnType],
   "Interpolation", Interpolation@Transpose@{ts, izTs},
   "TemporalData", TemporalData[{izTs}, {ts}]
   ]
  ]

Options[izTInterpolation] = {
   tPoints -> 200,
   ReturnType -> "Interpolation"
   };


(* 
stopPoints = 
  allPointsDS // Select[#numZPoints > 0 &] // Select[#stop > 0 &];
initialPoint = 
  stopPoints // Select[#numZPoints > 0 &] // RandomChoice // 
    Part[#, "initialPoint"] & // Normal;
sol = getPoint[initialPoint, ReturnSolution -> True, 
    PrecisionGoal -> \[Infinity]]["solution"];
overviewPlot[sol] 
*)


(* ::Subsection:: *)
(*Scattering*)


(* ::Subsubsection:: *)
(*Px*)


pxPairs[dataset_] := (
  initialPx = dataset[All, "initialPoint"] // Map[#["px"] &] // Normal;
  finalPx = dataset[All, "stopPoints"] // Flatten // Map[#["px"] &] // Normal;
  Thread[{initialPx, finalPx}]
  )

pxPairsPlot[pairs_List] := (
  labels = {"Initial px", "Final px"};
  plotOptions = {AspectRatio -> 1, PlotRange -> All, 
    AxesLabel -> labels};

  listPlot = ListPlot[pairs, plotOptions];
  denHist = DensityHistogram[pairs, FrameLabel -> labels, ChartLegends -> Automatic];

  {
   listPlot,
   Show[denHist, listPlot]
   }
  )

pxPairsPlot[dataset_Dataset] := (
  pxPairsPlot[pxPairs[dataset]]
  )



(* ::Subsubsection:: *)
(*Pitch angle*)


paAssoc[point_Association] := (
  If[KeyExistsQ[point, "x"],
   xTemp = point["x"]];
  If[KeyExistsQ[point, "kx"],
   xTemp = point["kx"]/k];
  
  pxTemp = point["px"];
  zTemp = point["z"];
  pzTemp = point["pz"];
  pa[xTemp, zTemp, pxTemp, pzTemp]
)

paPairs[dataset_Dataset] := (
  initialPa = dataset[All, "initialPoint"] // Map[paAssoc] // Normal;
  finalPa = dataset[All, "stopPoints"] // Flatten // Map[paAssoc] // Normal;
  Thread[{initialPa, finalPa}]
  )

paPairsPlot[pairs_List] := (
  labels = {"Initial pitch angle", "Final pitch angle"};
  plotOptions = {AspectRatio -> 1, PlotRange -> All, 
    AxesLabel -> labels};
  
  listPlot = ListPlot[pairs, plotOptions];
  denHist = DensityHistogram[pairs, FrameLabel -> labels, ChartLegends -> Automatic];

  {
   listPlot,
   Show[denHist, listPlot]
   }
  )

paPairsPlot[dataset_Dataset] := (
  paPairsPlot[paPairs[dataset]]
  )
