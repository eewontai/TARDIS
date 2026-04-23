server <- function(input, output, session) {

  shinyjs::useShinyjs()  # Initialize shinyjs - enables javascript - define actions

  volumes = shinyFiles::getVolumes()()  # detect hard drives on computer

  # open a folder selection window to choose raw data location
  shinyFiles::shinyDirChoose(input, "dir", roots = volumes)

  # display folder path as text to confirm
  output$selectedDir <- renderText({
    if (!is.null(input$dir)) {
      paste("Selected Directory:", shinyFiles::parseDirPath(volumes,input$dir))
    }
  })

  # same for output directory
  shinyFiles::shinyDirChoose(input, "dir_out", roots = volumes)

  output$selectedDirOut <- renderText({
    if (!is.null(input$dir_out)) {
      paste("Selected Directory:", shinyFiles::parseDirPath(volumes,input$dir_out))
    }
  })




  # Reactive values to store targets
  targets <- reactiveValues()  # creates a special list that listens for changes

  # Function to create target list
  observeEvent(input$create_target_list, {  # only runs when 'create target list' button is clicked

    req(input$target_file)  # requires that a file is uploaded, if not, stop

    # load the database using the input from user GUI
    targets$targets <- createTargetList(
      input_directory_targets = input$target_file$datapath,
      pos_pattern = input$pos_pattern,
      neg_pattern = input$neg_pattern,
      polarity = input$polarity,
      ion_column = input$ion_column,
      columns_of_interest = unlist(strsplit(input$columns_of_interest, ", "))
    )
    # update button label to 'Target list created' for feedback
    updateActionButton(inputId = "create_target_list", label = "Target list created")
    output$create_target_list_output <- renderPrint({  # prints preview of the raw data structure
      targets$targets
    })
  })

  # input box for mass limit only appears if user checks a box (range==TRUE)
  output$mass_low <- renderUI({  # renders reactive HTML
    if(input$range == TRUE){
      # create numericInput box
      numericInput("mass_low","Mass range - lower limit",value = 67)
    } else {
      # no input box
      NULL
    }
  })

  output$mass_high <- renderUI({
    if(input$range == TRUE){
      numericInput("mass_high","Mass range - upper limit",value = 1000)
    } else {
      NULL
    }
  })



  # Function to run tardis_peaks (triggers when "Run" button is clicked)
  observeEvent(input$run_tardis_peaks, {
    # freeze screen, show loading animation
    shinybusy::show_modal_spinner(
      spin = "double-bounce",
      color = "#003B6F",
      text = "Processing...",
      session = shiny::getDefaultReactiveDomain()  # find session (connection between user and server) automatically
    )

    # parsing batches
    # It takes a text input like "1, 10, 11, 20" and
    # turns it into pairs: [1, 10] and [11, 20]. This
    # tells the program which samples belong to which
    # experimental batch.
    batch <- as.numeric(unlist(strsplit(input$batch_positions, ",")))
    result <- list()
    for (i in seq(1, length(batch), by = 2)) {
      result[[length(result) + 1]] <- c(batch[i], batch[i + 1])
    }

    # launch the algorithm: pass GUI settings to the processing function
    tardis_output <- tardisPeaks(
      file_path = shinyFiles::parseDirPath(volumes,input$dir),
      dbData = targets$targets,
      ppm = as.numeric(input$ppm),
      rtdev = as.numeric(input$rtdev),
      mass_range = c(input$mass_low,input$mass_high),
      polarity = input$polarity,
      output_directory = paste0(shinyFiles::parseDirPath(volumes,input$dir_out),'/'),
      plots_samples = input$plot_samples,
      plots_QC = input$plot_QCs,
      diagnostic_plots = input$diagnostic_plots,
      batch_positions =  result,
      QC_pattern = input$QC_pattern,
      sample_pattern = input$sample_pattern,
      rt_alignment = input$rt_alignment,
      int_std_id = unlist(strsplit(input$int_std_id, ",")),
      screening_mode = input$screening_mode,
      smoothing = input$smoothing,
      max_int_filter = input$max_int_filter,
      num_cores = input$num_cores,
      rt_mode = input$rt_mode
    )

    # after it's done, the loading spinner disappears and the button label changes
    # uses session info to know which user's screen to unfreeze
    shinybusy::remove_modal_spinner(session = getDefaultReactiveDomain())

    # edited - THIS CODE NEEDS TO BE PLACED AFTER PROCESSING RUN!
    updateActionButton(inputId = "run_tardis_peaks", label = "Processing done! You may close this window.")
  })
}
