#' Provides a sample dataframe for informative heatmap visualization
#'
#' This function returns a predefined dataframe that represents data suitable for
#' creating an informative heatmap. The dataframe includes values for three groups
#' and their associated p-values, with row names indicating gene identifiers.
#'
#' @return A dataframe with sample data for an informative heatmap.
#' @examples
#' informative_heatmap_df <- get_informative_heatmap_df()
#' @export
get_informative_heatmap_df <- function() {
  informative_heatmap_df <- structure(list(group1_value = c(0.42, 0.13, -0.32, -0.42),
                                 group2_value = c(-0.36, 0.41, 0.21, 0.32),
                                 group3_value = c(-0.17, -0.15, 0.12, 0.11),
                                 pvalue1 = c(0.01, 0.12, 0.12, 0.02),
                                 pvalue2 = c(0.07, 0.03, 0.21, 0.08),
                                 pvalue3 = c(0.07, 0.3, 0.06, 0.21)),
                            class = "data.frame", row.names = c("gene1", "gene2", "gene3", "gene4"))

  return(informative_heatmap_df)
}

#' Provides a sample dataframe for clear scatterplot visualization
#'
#' This function returns a predefined dataframe that represents data suitable for
#' creating a clear scatterplot. The dataframe includes log2 fold change values, p-values,
#' and regulation direction for a set of genes across different tissues and time points.
#'
#' @return A dataframe with sample data for a clear scatterplot.
#' @examples
#' clear_scatterplot_df <- get_clear_scatterplot_df()
#' @export
get_clear_scatterplot_df <- function() {
  clear_scatterplot <- structure(list(

    ProbeName =
       c("A_51_P439092...33363", "A_52_P384419...194714",
         "A_51_P349213...579102", "A_52_P308681...1099", "A_51_P506593...23052",
         "A_52_P472397...686", "A_52_P382040...187763", "A_51_P253224...445534",
         "A_51_P212515...165620", "A_51_P496400...947", "A_52_P78966...578761",
         "A_52_P392663...26869", "A_52_P1188270...2240", "A_52_P628127...579140",
         "A_52_P237927...613531", "A_52_P510215...577721", "A_52_P215876...180505",
         "A_51_P480928...6341", "A_51_P361675...202605", "A_52_P762901...413042",
         "A_52_P494930...603091", "A_52_P546135...6876", "A_51_P185584...439265",
         "A_52_P706621...580357", "A_51_P500044...430379", "A_51_P225764...165332",
         "A_52_P522190...191320", "A_52_P308863...580401", "A_52_P306236...579621",
         "A_51_P250459...916", "A_52_P510583...20587", "A_52_P1197518...433587",
         "A_51_P406365...424048", "A_52_P1099932...579283", "A_51_P379385...577270",
         "A_52_P126146...165493", "A_51_P259245...165689", "A_52_P154311...29697",
         "A_52_P384822...8667", "A_52_P1005512...22540", "A_52_P476075...600667",
         "A_51_P490100...595446", "A_51_P492648...579657", "A_52_P643798...431268",
         "A_51_P304771...579957", "A_52_P168482...170414", "A_51_P380861...577408",
         "A_52_P126274...6392", "A_52_P330424...10301", "A_51_P467051...455",
         "A_52_P410383...832", "A_51_P426708...617071", "A_51_P263407...140",
         "A_51_P186703...180494", "A_51_P171288...40782", "A_51_P305896...580172",
         "A_52_P205772...186545", "A_51_P209122...16973", "A_52_P260010...596635",
         "A_52_P391110...590898", "A_52_P346256...188555", "A_51_P117924...5185",
         "A_51_P292757...170362", "A_52_P3825...1166", "A_51_P315634...579528",
         "A_51_P241769...577728", "A_52_P176245...583563", "A_52_P105599...31567",
         "A_52_P617542...607027", "A_52_P408690...579590", "A_51_P446645...880",
         "A_51_P511833...412412", "A_52_P14526...183027", "A_52_P655532...578307",
         "A_51_P166941...580226", "A_52_P400905...578526", "A_52_P338887...20964",
         "A_52_P715745...610467", "A_51_P239884...185349", "A_52_P114722...614872",
         "A_51_P200915...15175", "A_52_P458177...578041", "A_52_P16016...38229",
         "A_52_P558843...576474", "A_52_P309890...165937", "A_51_P106269...580075",
         "A_51_P130567...2367", "A_51_P381611...444951", "A_52_P268269...6838",
         "A_52_P209184...37939", "A_51_P515314...170102", "A_52_P389965...412355",
         "A_51_P242767...592637", "A_51_P241943...19460", "A_51_P269793...580531",
         "A_52_P325819...443708", "A_52_P145059...13984", "A_51_P256759...614276",
         "A_51_P492676...1264", "A_52_P548470...172211"),

    log2fc = c(-0.174723077,
         -0.247694685, 1.299996406, 1.747401972, -0.365981518, 1.758912918,
         0.3732152, -0.054890314, -1.644842352, 1.795242973, -1.281928252,
         0.313349559, -2.033641911, -1.083198709, 0.099148363, -0.986068049,
         -0.489574144, 0.936977674, -0.037236899, -0.804724058, -0.324581202,
         -1.091747125, 0.319718467, -1.48616533, 0.231022841, -1.303340755,
         -0.238381463, -0.82654247, -1.508721908, 1.943426417, -0.309243814,
         -0.255133801, 0.32144603, -0.845592489, 1.191034535, 1.691143015,
         1.475659486, 0.420332856, -1.135408962, -0.362487601, 0.386132095,
         0.360038532, -1.453722762, 0.236883933, -1.116479626, 0.978147341,
         -1.392903245, -0.746072838, -0.762224323, -2.618483941, -2.069264904,
         0.011088379, -1.541077788, 0.439340644, -0.007012264, 1.62679202,
         -0.357608144, -0.578184061, 0.384087067, -0.687779723, -0.575064065,
         -1.337148428, -0.641746266, 1.355403972, 1.960034466, 1.994427731,
         0.78709363, -0.146659005, -0.258031497, -1.062129633, 2.558087528,
         -1.256226278, 0.516280095, -1.55384002, 1.209427388, 1.846922243,
         0.727951654, 0.094892328, -0.239954966, 0.032723827, 0.89932695,
         -1.212533703, -0.046929769, -3.010833911, 1.911102145, 1.137606007,
         -1.673734877, 0.156257118, 1.009441921, -0.10195841, 0.936606515,
         0.787053035, 0.77800592, 0.390994634, 0.849282125, -0.070190259,
         -0.615548341, -0.036530906, -1.326249433, -0.822141261),

    p = c(0.791913617,
        0.764354557, 0.024774591, 0.045937321, 0.531668669, 0.035130042,
        0.605194496, 0.832965966, 0.033283207, 0.042537024, 0.020731403,
        0.6259764, 0.070225812, 0.025123353, 0.882410435, 0.010081838,
        0.429268519, 0.156622817, 0.931446705, 0.055514113, 0.567513865,
        0.167967381, 0.688340289, 0.041010476, 0.477636216, 0.022604889,
        0.68916475, 0.041586358, 0.031427603, 0.042054287, 0.471875219,
        0.555014385, 0.333164369, 0.026671507, 0.006406498, 0.028140733,
        0.03568997, 0.696597179, 0.205178177, 0.518412894, 0.498642861,
        0.361338785, 0.031876774, 0.49867289, 0.035812217, 0.181231559,
        0.007588725, 0.157767285, 0.240872913, 0.027771825, 0.039526205,
        0.984779795, 0.012241003, 0.429011546, 0.989007891, 0.038545102,
        0.576265507, 0.388197423, 0.391336552, 0.248992465, 0.624051935,
        0.132688898, 0.179918758, 0.047339372, 0.030027321, 0.010115655,
        0.097069701, 0.744353687, 0.684869115, 0.030929349, 0.041093571,
        0.033776928, 0.490965389, 0.01588564, 0.0391728, 0.018244263,
        0.480288646, 0.788308656, 0.547675675, 0.922206072, 0.34703949,
        0.013100026, 0.918609426, 0.000125104, 0.044525656, 0.037221981,
        0.072952467, 0.819136534, 0.167288906, 0.910951548, 0.172494378,
        0.031441156, 0.291533955, 0.444958688, 0.043649914, 0.790183777,
        0.320995909, 0.905079639, 0.049593468, 0.227596452),

    regulation = c("down",
        "down", "up", "up", "down", "up", "up", "down", "down", "up",
        "down", "up", "down", "down", "up", "down", "down", "up", "down",
        "down", "down", "down", "up", "down", "up", "down", "down", "down",
        "down", "up", "down", "down", "up", "down", "up", "up", "up",
        "up", "down", "down", "up", "up", "down", "up", "down", "up",
        "down", "down", "down", "down", "down", "up", "down", "up", "down",
        "up", "down", "down", "up", "down", "down", "down", "down", "up",
        "up", "up", "up", "down", "down", "down", "up", "down", "up",
        "down", "up", "up", "up", "up", "down", "up", "up", "down", "down",
        "down", "up", "up", "down", "up", "up", "down", "up", "up", "up",
        "up", "up", "down", "down", "down", "down", "down"),

    organ = c("AY",
       "AY", "HC", "AY", "AY", "AY", "AY", "HC", "AY", "AY", "HC", "AY",
       "AY", "HC", "HC", "HC", "AY", "AY", "AY", "HC", "HC", "AY", "HC",
       "HC", "HC", "AY", "AY", "HC", "HC", "AY", "AY", "HC", "HC", "HC",
       "HC", "AY", "AY", "AY", "AY", "AY", "HC", "HC", "HC", "HC", "HC",
       "AY", "HC", "AY", "AY", "AY", "AY", "HC", "AY", "AY", "AY", "HC",
       "AY", "AY", "HC", "HC", "AY", "AY", "AY", "AY", "HC", "HC", "HC",
       "AY", "HC", "HC", "AY", "HC", "AY", "HC", "HC", "HC", "AY", "HC",
       "AY", "HC", "AY", "HC", "AY", "HC", "AY", "HC", "AY", "HC", "AY",
       "AY", "AY", "HC", "HC", "AY", "HC", "HC", "AY", "HC", "AY", "AY"),

    timePoint =
      structure(
        c(1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 2L,
         1L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 1L, 2L, 1L,
         2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 2L,
         2L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 2L,
         1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 2L,
         2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 1L,
         1L, 2L, 1L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L),

        levels = c("T10R1", "T5R1"),

        class = "factor"),

    reg_time_org = c("down.AY.T10R1",
        "down.AY.T5R1", "up.HC.T5R1", "up.AY.T10R1", "down.AY.T10R1",
        "up.AY.T10R1", "up.AY.T5R1", "down.HC.T10R1", "down.AY.T5R1",
        "up.AY.T10R1", "down.HC.T5R1", "up.AY.T10R1", "down.AY.T10R1",
        "down.HC.T5R1", "up.HC.T5R1", "down.HC.T5R1", "down.AY.T5R1",
        "up.AY.T10R1", "down.AY.T5R1", "down.HC.T10R1", "down.HC.T5R1",
        "down.AY.T10R1", "up.HC.T10R1", "down.HC.T5R1", "up.HC.T10R1",
        "down.AY.T5R1", "down.AY.T5R1", "down.HC.T5R1", "down.HC.T5R1",
        "up.AY.T10R1", "down.AY.T10R1", "down.HC.T10R1", "up.HC.T10R1",
        "down.HC.T5R1", "up.HC.T5R1", "up.AY.T5R1", "up.AY.T5R1", "up.AY.T10R1",
        "down.AY.T10R1", "down.AY.T10R1", "up.HC.T5R1", "up.HC.T5R1",
        "down.HC.T5R1", "up.HC.T10R1", "down.HC.T5R1", "up.AY.T5R1",
        "down.HC.T5R1", "down.AY.T10R1", "down.AY.T10R1", "down.AY.T10R1",
        "down.AY.T10R1", "up.HC.T5R1", "down.AY.T10R1", "up.AY.T5R1",
        "down.AY.T10R1", "up.HC.T5R1", "down.AY.T5R1", "down.AY.T10R1",
        "up.HC.T5R1", "down.HC.T5R1", "down.AY.T5R1", "down.AY.T10R1",
        "down.AY.T5R1", "up.AY.T10R1", "up.HC.T5R1", "up.HC.T5R1", "up.HC.T5R1",
        "down.AY.T10R1", "down.HC.T5R1", "down.HC.T5R1", "up.AY.T10R1",
        "down.HC.T10R1", "up.AY.T5R1", "down.HC.T5R1", "up.HC.T5R1",
        "up.HC.T5R1", "up.AY.T10R1", "up.HC.T5R1", "down.AY.T5R1", "up.HC.T5R1",
        "up.AY.T10R1", "down.HC.T5R1", "down.AY.T10R1", "down.HC.T5R1",
        "up.AY.T5R1", "up.HC.T5R1", "down.AY.T10R1", "up.HC.T10R1", "up.AY.T10R1",
        "down.AY.T10R1", "up.AY.T5R1", "up.HC.T10R1", "up.HC.T5R1", "up.AY.T10R1",
        "up.HC.T5R1", "down.HC.T10R1", "down.AY.T10R1", "down.HC.T5R1",
        "down.AY.T10R1", "down.AY.T5R1"),

    negLog10p = c(0.101322189251776,
        0.116705140861921, 1.60599350648585, 1.33783433584631, 0.274358931713898,
        1.45432133098346, 0.218105030364148, 0.0793727429873967, 1.47777483396483,
        1.37123289773345, 1.68337130598429, 0.203442039861418, 1.15350323048052,
        1.59992239948019, 0.0543293647608365, 1.99646028530209, 0.367270959938686,
        0.805144969127632, 0.0308419892772071, 1.25559659491582, 0.246023523704117,
        0.774775049388113, 0.162196810166872, 1.38710518991823, 0.320902750947676,
        1.64579762120239, 0.161676944329053, 1.38104911208905, 1.50268874159047,
        1.37618972577465, 0.326172829468126, 0.255695760579735, 0.477341451365263,
        1.57395244500401, 2.19337930509684, 1.55066459440272, 1.4474538171001,
        0.157018288630602, 0.687868833192218, 0.285324205201873, 0.302210394312456,
        0.442085420088633, 1.49652563661376, 0.30218424125313, 1.44596879261381,
        0.741766173623127, 2.11983118483946, 0.801983048128944, 0.618212035186559,
        1.55639558004856, 1.40311488124715, 0.00666087052558755, 1.91218299562845,
        0.367531019478186, 0.00480024328241102, 1.41403080073932, 0.239377374810247,
        0.410947351947642, 0.407449585806132, 0.603813795307619, 0.204779265853875,
        0.877165412772641, 0.744923555547126, 1.32477750768162, 1.52248341328278,
        1.99500599092819, 1.01291630811411, 0.128220656102862, 0.16439241838276,
        1.50962922091771, 1.38622611723878, 1.47137985190287, 0.308949122736756,
        1.79899528364644, 1.40701538476848, 1.738873675818, 0.318497679939908,
        0.103303704653155, 0.261476547820299, 0.0351720226203927, 0.459621103557027,
        1.88272784238658, 0.0368691025292607, 3.90272880421421, 1.35138967332274,
        1.42920051689912, 1.13696001718559, 0.0866437038308143, 0.77653285893442,
        0.0405047218149637, 0.763225055045001, 1.50250149459218, 0.535310855620789,
        0.351680309033674, 1.3600166076139, 0.102271890918435, 0.493500502517393,
        0.04331320504172, 1.30457552105753, 0.642834512438795),

    color_flag = c(0,
         0, 1, 1, 0, 1, 0, 0, -1, 1, -1, 0, 0, -1, 0, -1, 0, 0, 0, 0,
         0, 0, 0, -1, 0, -1, 0, -1, -1, 1, 0, 0, 0, -1, 1, 1, 1, 0, 0,
         0, 0, 0, -1, 0, -1, 0, -1, 0, 0, -1, -1, 0, -1, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, -1, 1, -1, 0, -1, 1, 1, 0, 0,
         0, 0, 0, -1, 0, -1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
         -1, 0)),

    row.names = c(NA, -100L),

    class = "data.frame")


  return(clear_scatterplot)
}

#' Provides a sample dataframe for multifeature grid visualization
#'
#' This function returns a predefined dataframe that represents data suitable for
#' creating a multifeature grid visualization. The dataframe includes signaling pathways,
#' activation z-scores, and associated p-values for different tissues and time points.
#'
#' @return A dataframe with sample data for a multifeature grid visualization.
#' @examples
#' multifeature_grid_df <- get_multifeature_grid_df()
#' @export
get_multifeature_grid_df <- function() {
  multifeature_grid <- structure(list(signaling = c("pathway1", "pathway2", "pathway3",
                                                    "pathway4", "pathway1", "pathway2", "pathway3", "pathway4", "pathway1",
                                                    "pathway2", "pathway3", "pathway4", "pathway1", "pathway2", "pathway3",
                                                    "pathway4", "pathway1", "pathway2", "pathway3", "pathway4", "pathway1",
                                                    "pathway2", "pathway3", "pathway4", "pathway1", "pathway2", "pathway3",
                                                    "pathway4", "pathway1", "pathway2", "pathway3", "pathway4", "pathway1",
                                                    "pathway1", "pathway2", "pathway2", "pathway3", "pathway3", "pathway4",
                                                    "pathway4"),
                                      Activation_z_score = c(-1.37, -0.34, -2.227, -1.642,
                                              -2.15, 1.352, 0.215, 0.013, 1.37, -2.193, -0.738, 1.433, 1.37,
                                              1.268, -0.563, 0.562, 0.033, -1.974, -1.041, 0.247, 1.37, -0.359,
                                              -0.6, -0.285, 1.37, -0.905, -1.516, -1.606, 0.813, 0.897, 0.766,
                                              0.937, NA, NA, NA, NA, NA, NA, NA, NA),
                                      timePoint = c("  acute  stress",
                                              "  acute  stress", "  acute  stress", "  acute  stress", "  acute  stress",
                                              "  acute  stress", "  acute  stress", "  acute  stress", "  acute recovery",
                                              "  acute recovery", "  acute recovery", "  acute recovery", "  acute recovery",
                                              "  acute recovery", "  acute recovery", "  acute recovery", " chronic  stress",
                                              " chronic  stress", " chronic  stress", " chronic  stress", " chronic  stress",
                                              " chronic  stress", " chronic  stress", " chronic  stress", " chronic recovery",
                                              " chronic recovery", " chronic recovery", " chronic recovery",
                                              " chronic recovery", " chronic recovery", " chronic recovery",
                                              " chronic recovery", "DEG", "DEG", "DEG", "DEG", "DEG", "DEG",
                                              "DEG", "DEG"),
                                      tissue = c("tissue1", "tissue1", "tissue1", "tissue1",
                                              "tissue2", "tissue2", "tissue2", "tissue2", "tissue1", "tissue1",
                                              "tissue1", "tissue1", "tissue2", "tissue2", "tissue2", "tissue2",
                                              "tissue1", "tissue1", "tissue1", "tissue1", "tissue2", "tissue2",
                                              "tissue2", "tissue2", "tissue1", "tissue1", "tissue1", "tissue1",
                                              "tissue2", "tissue2", "tissue2", "tissue2", "genes", "genes",
                                              "genes", "genes", "genes", "genes", "genes", "genes"),
                                      p = c(NA,
                                               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 4.49e-42,
                                               1.15e-42, 1.64e-20, 1.42e-21, 6.03e-21, 1.73e-22, 7.95e-30, 9.14e-30
                                              ),
                                      number_of_genes = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                               NA, NA, NA, NA, NA, NA, 40L, 40L, 26L, 27L, 23L, 24L, 40L, 40L
                                          )),
                                 class = "data.frame", row.names = c(NA, -40L))

  return(multifeature_grid)
}
