ui <- bslib::page_navbar(

  shinyjs::useShinyjs(),  # Initialize shinyjs
  title = div(img(src="tardis_new.png",height=100,width =100)),
  bg = "#003B6F",
  inverse = TRUE,  # background is dark, text is light color
  # 'create target list' tab - select/input options in boxes
  bslib::nav_panel(title = "Create target list",
            p(
              column(width = 6,
                     fileInput("target_file", "Table with target compounds (.csv or .xlsx):", accept = c(".xlsx",".csv")),
                     textInput("pos_pattern", "Positive pattern:", value = "+"),
                     textInput("neg_pattern", "Negative pattern:", value = "-"),
                     selectInput("polarity", "Polarity:", choices = c("positive","negative")),
                     textInput("ion_column", "Name of ion column:", value = "ion"),
                     textInput("columns_of_interest", "Names of columns of interest (comma-separated):", value = "id, name, mz, rt"),
                     actionButton("create_target_list", "Create target list")
              )
            )),
  # 'targeted peak detection' tab - select/input options in boxes
  bslib::nav_panel(title = "Targeted peak detection",
            bslib::page_fillable(
              bslib::layout_columns(
                bslib::card(
                  shinyFiles::shinyDirButton("dir", "Enter data folder path:",title ="directory "),
                  textOutput("selectedDir"),
                  shinyFiles::shinyDirButton("dir_out", "Enter output folder path:",title ="directory out"),
                  textOutput("selectedDirOut")
                ),
                bslib::card(
                  numericInput("ppm", "PPM:", value = 5),
                  #selectInput("mode", "Mode:", choices = c("metabolomics","lipidomics")),
                  checkboxInput("range", "Custom mass range: needed for data with >1 scan window!", value = TRUE),
                  uiOutput("mass_low"),
                  uiOutput("mass_high"),
                  textInput("rtdev", "RT deviation:", value = "18"),
                  textInput("int_std_id", "Internal standard IDs (comma-separated):", value = "331,1578,1576,1583,1577"),
                  textInput("batch_positions", "Batch postions (comma-separated):", value = "1,30"),
                  textInput("sample_pattern", "Sample pattern:", value = ""),
                  textInput("QC_pattern", "QC pattern:", value = "QC")
                ),
                bslib::card(
                  checkboxInput("diagnostic_plots", "Generate diagnostic plots", value = TRUE)
                  checkboxInput("plot_samples", "Generate sample plots", value = FALSE),
                  checkboxInput("plot_QCs", "Generate QC plots", value = FALSE),
                  checkboxInput("screening_mode", "Screening mode", value = FALSE),
                  checkboxInput("rt_alignment", "RT alignment", value = TRUE),
                  checkboxInput("smoothing", "Smoothing", value = TRUE),
                  numericInput("max_int_filter", "Peak intensity filter:", value = NULL),
                  numericInput("num_cores", "Number of Cores to use for Parallelization:", value = 1),
                  selectInput("rt_mode", "Retention time alignment mode:", choices = c("mean","linear_interpolation")),
                  numericInput("pval_cutoff", "P-value cutoff for diptest:", value = 0.05)
                ),
                bslib::card(
                  actionButton("run_tardis_peaks", "Run TARDIS"),


                )
              )
            ))


)
