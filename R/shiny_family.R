#' @title  plotSurfaceInteractiveSingleCell
#' @author Dieter Henrik Heiland
#' @description plotSurfaceInteractiveSingleCell
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotSurfaceInteractiveSingleCell <- function(object){
    SPATA2::check_object(object)
    surface_plots <- shiny::runApp(shiny::shinyApp(ui = function() {
      shinydashboard::dashboardPage(shinydashboard::dashboardHeader(title = "Surface Plots"), 
                                    shinydashboard::dashboardSidebar(collapsed = TRUE, 
                                                                     shinydashboard::sidebarMenu(shinydashboard::menuItem(text = "Surface Plots", 
                                                                                                                          tabName = "surface_plots", selected = TRUE))), 
                                    shinydashboard::dashboardBody(shinybusy::add_busy_spinner(spin = "cube-grid", 
                                                                                              margins = c(0, 10), color = "red"), shinydashboard::tabItems(tab_surface_plots_return())))
    }, server = function(input, output, session) {
      output$saved_plots <- shiny::renderUI({
        saved_plots <- base::names(plot_list())
        shiny::validate(shiny::need(base::length(saved_plots) != 
                                      0, message = "No plots have been saved yet."))
        shinyWidgets::checkboxGroupButtons(inputId = "saved_plots", 
                                           label = "Choose plots to export", choices = saved_plots, 
                                           selected = saved_plots, checkIcon = list(yes = icon("ok", 
                                                                                               lib = "glyphicon")))
      })
      plot_list <- shiny::reactiveVal(value = list())
      plot_df <- shiny::reactiveVal(value = data.frame())
      module_return <- moduleSurfacePlotServer(id = "isp", 
                                               object = object, final_plot = shiny::reactive(module_return()$assembled_plot()), 
                                               reactive_object = shiny::reactive(object))
      final_plot <- shiny::reactive({
        module_return()$assembled_plot()
      })
      oe <- shiny::observeEvent(input$save_plot, {
        plot_list <- plot_list()
        if (input$plot_name %in% base::names(plot_list) | 
            input$plot_name == "") {
          shiny::showNotification(ui = "Plot name is already taken or invalid.", 
                                  type = "error")
        }
        else {
          plot_list[[input$plot_name]] <- final_plot()
          plot_list(plot_list)
          shiny::showNotification(ui = "Plot has been saved.", 
                                  type = "message")
        }
      })
      oe <- shiny::observeEvent(input$return_plot, {
        plot_list <- plot_list()
        shiny::stopApp(returnValue = plot_list[base::names(plot_list) %in% 
                                                 input$saved_plots])
      })
      output$surface_variable <- shiny::renderPlot({
        plot_df <- module_return()$smoothed_df()
        var_name <- base::colnames(plot_df)[5]
        if (base::is.numeric(dplyr::pull(plot_df, var_name))) {
          plot_type <- input$surface_variable_plot_type
          if (plot_type == "violin") {
            add_on <- ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                                     axis.ticks.x = ggplot2::element_blank())
          }
          else {
            add_on <- list()
          }
          plotDistribution2(df = plot_df, plot_type = plot_type, 
                            binwidth = 0.05, verbose = FALSE) + add_on
        }
        else {
          ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = .data[[var_name]])) + 
            ggplot2::geom_bar(mapping = ggplot2::aes(fill = .data[[var_name]]), 
                              color = "black") + ggplot2::theme_classic() + 
            ggplot2::theme(legend.position = "none") + 
            confuns::scale_color_add_on(aes = "fill", variable = "discrete", 
                                        clrp = module_return()$current_setting()$pt_clrp) + 
            ggplot2::labs(y = "Count")
        }
      })
    }))
    return(surface_plots)
  
}


#' @title Surface plot tab - return
#' @details To use within shinydashboard::tab_items()
#' @note Tab for the output returning application
#'
#' @export

tab_surface_plots_return <- function(){shinydashboard::tabItem(tabName = "surface_plots",
                                                               
                                                               shiny::fluidRow(
                                                                 shiny::column(width = 7, align = "center",
                                                                               moduleSurfacePlotUI(id = "isp"),
                                                                               
                                                                 ),
                                                                 
                                                                 shiny::column(width = 5, align = "center",
                                                                               shiny::wellPanel(
                                                                                 shiny::fluidRow(
                                                                                   shiny::column(width = 12,
                                                                                                 shiny::plotOutput("surface_variable"),
                                                                                                 shiny::HTML("<br>"),
                                                                                                 shinyWidgets::radioGroupButtons(
                                                                                                   inputId = "surface_variable_plot_type",
                                                                                                   label = NULL,
                                                                                                   selected = "density",
                                                                                                   choices = c("Densityplot" = "density",
                                                                                                               "Histogram" = "histogram",
                                                                                                               "Violinplot" = "violin")
                                                                                                 )
                                                                                   )
                                                                                 )
                                                                               )
                                                                 )
                                                                 
                                                               ),
                                                               shiny::fluidRow(
                                                                 shiny::column(width = 4, align = "center",
                                                                               shiny::textInput("plot_name", label = NULL, value = "", placeholder = "Plot name"),
                                                                               shiny::actionButton("save_plot", label = "Save Plot"),
                                                                               shiny::actionButton("return_plot", label = "Return Plots")
                                                                 ),
                                                                 shiny::column(width = 1, align = "center",
                                                                               shiny::uiOutput("saved_plots")
                                                                 )
                                                               )
)}


#' @title Surface plot tab - classic
#' @details To use within shinydashboard::tab_items()
#' @note Tab for the big application#' 
#' @export

tab_surface_plots_app <- function(){shinydashboard::tabItem(tabName = "surface_plots",
                                                            
                                                            shiny::fluidRow(
                                                              shiny::column(width = 7, align = "center",
                                                                            moduleSurfacePlotUI(id = "isp"),
                                                                            
                                                              ),
                                                              
                                                              shiny::column(width = 5, align = "center",
                                                                            shiny::wellPanel(
                                                                              shiny::fluidRow(
                                                                                shiny::column(width = 12,
                                                                                              shiny::plotOutput("surface_variable"),
                                                                                              shiny::HTML("<br>"),
                                                                                              shinyWidgets::radioGroupButtons(
                                                                                                inputId = "surface_variable_plot_type",
                                                                                                label = NULL,
                                                                                                selected = "density",
                                                                                                choices = c("Densityplot" = "density",
                                                                                                            "Histogram" = "histogram",
                                                                                                            "Violinplot" = "violin")
                                                                                              )
                                                                                )
                                                                              )
                                                                            )
                                                              )
                                                              
                                                            )
)}

#' @title UI of the surface plot module
#'
#' @param id The namespace id.
#'
#' @export

moduleSurfacePlotUI <- function(id){
  
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::column(width = 12,
                  shiny::wellPanel(
                    shiny::fluidRow(
                      shiny::column(width = 4,
                                    shiny::fluidRow(
                                      shiny::column(width = 6, shiny::uiOutput(ns("sample_opts"))),
                                      shiny::column(width = 6, shiny::uiOutput(ns("aes_clr_opts")))
                                    ),
                                    shiny::fluidRow(
                                      shiny::column(width = 12,
                                                    shiny::uiOutput(ns("aes_clr_opts_detailed")),
                                                    shiny::conditionalPanel(
                                                      condition = "input.aes_clr_opts == 'gene_sets'", ns = ns,
                                                      shinyWidgets::pickerInput(ns("method_gs"),
                                                                                label = "Gene-set method:",
                                                                                choices = c("Mean" = "mean",
                                                                                            "Gene Set Variation Analysis" = "gsva",
                                                                                            "Gene Set Enrichment Analysis" = "ssgsea",
                                                                                            "Z-Score" = "zscore",
                                                                                            "Plage" = "plage" )))),
                                    ),
                                    shiny::fluidRow(
                                      shiny::column(width = 6, shiny::uiOutput(ns("pt_clrsp"))),
                                      shiny::column(width = 6, shiny::uiOutput(ns("pt_clrp")))
                                    ),
                                    shiny::fluidRow(
                                      shiny::column(width = 6,
                                                    shiny::sliderInput(ns("pt_size"), label = "Size of points:", min = 0.1, max = 10, step = 0.05, value = 2.75),
                                                    shiny::sliderInput(ns("pt_alpha"), label = "Transparency of points:", min = 0.01, max = 0.99, step = 0.01, value = 0.15),
                                                    shiny::uiOutput(ns("pt_smooth"))
                                      ),
                                      shiny::column(width = 6,
                                                    shiny::uiOutput(ns("scale_color_min")),
                                                    shiny::uiOutput(ns("scale_color_mid")),
                                                    shiny::uiOutput(ns("scale_color_max"))
                                      )
                                    ),
                                    shiny::HTML("<br>")
                      ),
                      shiny::column(width = 8,
                                    shiny::plotOutput(ns("surface_plot"), dblclick = ns("surface_plot_dblclick")),
                                    shiny::HTML("<br>"),
                                    shiny::fluidRow(
                                      shiny::column(width = 4,
                                                    shiny::actionButton(ns("update_plot"), label = "Plot & Update")),
                                      shiny::column(width = 8,
                                                    shinyWidgets::checkboxGroupButtons(inputId = ns("display_add_ons"),
                                                                                       label = NULL,
                                                                                       selected = "legend",
                                                                                       choices = c("Legend" = "legend",
                                                                                                   "Image" = "image",
                                                                                                   "Title" = "title",
                                                                                                   "Coordinates" = "coords",
                                                                                                   "Segmentation" = "segmentation"),
                                                                                       direction = "horizontal",
                                                                                       justified = FALSE,
                                                                                       individual = FALSE)
                                      )
                                    )
                      )
                    )
                  )
    )
  )
  
}


#' @title Server of the surface plot module
#'
#' @param id  The namespace id.
#' @param object A valid spata-object.
#' @param final_plot The final plot that is to be displayed. (See details.).
#' @param reactive_object A valid (reactive) spata-object.
#'
#' @return A reactive list with several slots:
#'  \enumerate{
#'   \item $assembled_plot() The surface plot as a ggplot-object.
#'   \item $dblclick() A list containing information regarding the double clicked position in the plot.
#'   \item $current_setting() A list with information about the settings of \code{assembled_plot} (e.g. sample, color_to, smooth, smoothing_span ...)}
#'
#' @details The argument \code{final_plot} takes a ggplot object as input which is going to be displayed as the final plot. This allows to
#' adjust the output of \code{$assembled_plot()} outside of the module. If no further adjustment is needed determine \code{final_plot} as:
#' \code{shiny::reactive(*module_return_variable*()$assembled_plot())}
#' @export

moduleSurfacePlotServer <- function(id,
                                    object,
                                    final_plot,
                                    reactive_object,
                                    highlighted = shiny::reactive( FALSE )){
  
  shiny::moduleServer(
    id = id,
    module = function(input,
                      output,
                      session){
      
      # Reactive values -----------------------------------------------------------
      
      return_plot <- shiny::reactiveVal(list())
      
      current <- shiny::reactiveValues(
        
        sample = getSampleNames(object)[1],
        color_code = "gene_sets",
        gene_set = base::character(1),
        method_gs = base::character(1),
        genes = base::character(1),
        feature = base::character(1),
        pt_size = base::numeric(1),
        pt_clrp = base::character(1),
        pt_clrsp = base::character(1),
        smooth = base::logical(1),
        span = base::numeric()
        
      )
      
      reset_select_gene_sets <- shiny::reactiveVal(value = 0)
      reset_select_genes <- shiny::reactiveVal(value = 0)
      
      all_features <- getFeatureNames(object) %>% base::unname()
      all_gene_sets <- getGeneSets(object = object)
      all_genes <- getGenes(object = object, in_sample = "all")
      
      smooth_values <- base::seq(0.01, 0.25, by = 0.01) %>%
        base::round(digits = 3) %>%
        base::unique()
      
      all_values <- c(0, smooth_values)
      
      # -----
      
      # Render UIs and Outputs --------------------------------------------------
      
      # update transparency
      
      shiny::observeEvent(eventExpr = highlighted(), {
        
        if(base::isTRUE(highlighted())){
          
          shiny::updateSliderInput(session,
                                   inputId = "pt_alpha",
                                   label = "Transparency of points",
                                   min = 0.01,
                                   max = 0.99,
                                   step = 0.01,
                                   value = 0.75)
          
        } else if(base::isFALSE(highlighted())){
          
          shiny::updateSliderInput(session,
                                   inputId = "pt_alpha",
                                   label = "Transparency of points",
                                   min = 0.01,
                                   max = 0.99,
                                   step = 0.01,
                                   value = 0.15)
          
        }
        
      })
      
      # Main select input -------------------------------------------------------
      
      output$sample_opts <- shiny::renderUI({
        
        ns <- session$ns
        
        shinyWidgets::pickerInput(ns("sample_opts"),
                                  label = "Choose sample:",
                                  choices = getSampleNames(object),
                                  selected = getSampleNames(object)[1])
        
      })
      
      output$aes_clr_opts <- shiny::renderUI({
        
        ns <- session$ns
        
        shinyWidgets::pickerInput(ns("aes_clr_opts"),
                                  label = "Color by:",
                                  choices = c("Gene set" = "gene_sets",
                                              "Genes" = "genes",
                                              "Feature" = "feature"),
                                  selected = "feature")
        
      })
      
      select_gene_sets <- shiny::eventReactive(reset_select_gene_sets(),{
        
        ns <- session$ns
        
        shinyWidgets::pickerInput(inputId = ns("aes_clr_opts_detailed"),
                                  label = "Choose gene-set:",
                                  choices = all_gene_sets,
                                  selected = all_gene_sets[1],
                                  options = list(`live-search` = TRUE),
                                  multiple = F)
        
      })
      
      select_genes <- shiny::eventReactive(reset_select_genes(),{
        
        ns <- session$ns
        
        shiny::tagList(
          shinyWidgets::pickerInput(inputId = ns("aes_clr_opts_detailed"),
                                    label = "Choose gene(s):",
                                    choices = all_genes,
                                    selected = all_genes[1],
                                    options = shinyWidgets::pickerOptions(
                                      liveSearch = TRUE,
                                      actionsBox = TRUE),
                                    multiple = TRUE),
          shiny::checkboxInput(ns("reset_select_genes"),
                               label = "Automatic reset",
                               value = FALSE))
        
      })
      
      select_features <- shiny::reactive({
        
        ns <- session$ns
        
        shinyWidgets::pickerInput(inputId = ns("aes_clr_opts_detailed"),
                                  label = "Choose feature:",
                                  choices = all_features[all_features != "sample"],
                                  options = shinyWidgets::pickerOptions(
                                    liveSearch = TRUE,
                                    actionsBox = TRUE),
                                  multiple = F)
        
      })
      
      output$aes_clr_opts_detailed <- shiny::renderUI({
        
        shiny::req(input$aes_clr_opts)
        
        if(input$aes_clr_opts == "gene_sets"){
          
          base::return(select_gene_sets())
          
        } else if(input$aes_clr_opts == "genes"){
          
          base::return(select_genes())
          
        } else if(input$aes_clr_opts == "feature"){
          
          base::return(select_features())
          
        }
        
      })
      
      # -----
      
      # Color select input ------------------------------------------------------
      
      output$pt_clrsp <- shiny::renderUI({
        
        ns <- session$ns
        
        shinyWidgets::pickerInput(ns("pt_clrsp"),
                                  label = "Color spectrum:",
                                  choices = validColorSpectra(),
                                  options = list(
                                    `live-search` = TRUE
                                  ),
                                  multiple = FALSE,
                                  selected = "inferno")
        
      })
      
      output$pt_clrp <- shiny::renderUI({
        
        ns <- session$ns
        
        choices = c(
          "MILO Research Group" = "milo",
          "Journal of Oncology" = "jco",
          "Nature Publishing Group" = "npg",
          "American Association for the Advancement" = "aaas",
          "New England Journal of Medicine" = "nejm",
          "Lancet Oncology" = "lo",
          "The Journal of the American Medical Association" = "jama",
          "University of Chicago" = "uc")
        
        shinyWidgets::pickerInput(ns("pt_clrp"),###!
                                  choices = validColorPalettes(),
                                  label = "Color palette:",
                                  multiple = FALSE,
                                  choicesOpt = list(
                                    #subtext = stringr::str_c("colors: ", c(20, base::rep(10,7))),
                                    `dropdown-align-center` = TRUE
                                  ),
                                  selected = "milo")
      })
      
      # -----
      
      
      
      # Plot tweaking slider inputs ---------------------------------------------
      
      output$scale_color_min <- shiny::renderUI({
        
        shiny::validate(
          shiny::need(base::is.numeric(color_variable()),
                      message = "Need numeric color-feature to scale minimum.",
                      label = "Color scale minimum")
        )
        
        ns <- session$ns
        
        shiny::sliderInput(ns("scale_color_min"),
                           label = "Color scale minimum:",
                           min = color_min(),
                           max = color_max(),
                           value = color_min(),
                           step = 0.01)
        
      })
      
      output$scale_color_max <- shiny::renderUI({
        
        shiny::validate(
          shiny::need(expr = base::is.numeric(color_variable()),
                      message = "Need numeric color-feature to scale maximum.",
                      label = "Color scale maximum:")
        )
        
        ns <- session$ns
        
        shiny::sliderInput(ns("scale_color_max"),
                           label = "Color scale maximum:",
                           min = color_min(),
                           max = color_max(),
                           value = color_max(),
                           step = 0.01)
        
      })
      
      output$scale_color_mid <- shiny::renderUI({
        
        shiny::req(base::is.numeric(color_variable()))
        
        ns <- session$ns
        
        shiny::sliderInput(ns("scale_color_mid"),
                           label = "Color scale mid:",
                           min = color_min() * 1.1,
                           max = color_max() * 0.9,
                           value = color_median(),
                           step = 0.01)
        
      })
      
      output$pt_smooth <- shiny::renderUI({
        
        ns <- session$ns
        
        shinyWidgets::sliderTextInput(
          inputId = ns("pt_smooth"),
          label = "Spatial smoothing:",
          choices = all_values,
          grid = TRUE,
          selected = 0
        )
        
      })
      
      # -----
      
      
      
      # Plot assembling ---------------------------------------------------------
      
      output$surface_plot <- shiny::renderPlot({
        
        shiny::req(final_plot())
        
        final_plot()
        
      })
      
      # -----
      
      # Plot add-ons ------------------------------------------------------------
      
      #----- Image add-on -----#
      
      image_add_on <- shiny::reactive({
        
        ## set up background
        if("image" %in% input$display_add_ons){
          
          ## extract image info
          img_info <-
            getImage(object, of_sample = current$sample) %>%
            grDevices::as.raster() %>%
            magick::image_read() %>%
            magick::image_info()
          
          st_image <-
            grDevices::as.raster(getImage(object, of_sample = current$sample)) %>%
            magick::image_read() %>% magick::image_flip()
          
          image_add_on <-
            ggplot2::annotation_raster(raster = st_image,
                                       xmin = 0, ymin = 0,
                                       xmax = img_info$width,
                                       ymax = img_info$height)
          
          
        } else {
          
          image_add_on <- NULL
          
        }
        
        
      })
      
      #----- Geom point add-on -----#
      
      # sample coordinates
      sample_coords <- shiny::reactive({
        
        sample_coords <-
          getCoordsDf(object = object, of_sample = current$sample)
        
        base::return(sample_coords)
        
      })
      
      # rna_assay
      rna_assay <- shiny::reactive({
        
        rna_assay <-
          getExpressionMatrix(object = object, of_sample = current$sample)
        
        base::return(rna_assay)
        
      })
      
      # gene_vls
      gene_vls <- shiny::reactive({
        
        genes <- current$genes
        
        # compute mean if neccessary
        if(base::length(genes) > 1){
          
          rna_assay <- base::colMeans(rna_assay()[genes,])
          
        } else {
          
          rna_assay <- rna_assay()[genes,]
          
        }
        
        
        # convert to data frame
        gene_vls <-
          rna_assay %>%
          as.data.frame() %>%
          magrittr::set_colnames(value = "expr_score") %>%
          tibble::rownames_to_column(var = "barcodes")
        
        base::return(gene_vls)
        
      })
      
      # geneset_vls
      geneset_vls <- shiny::reactive({
        
        shiny::req(current$gene_set)
        
        gene_set_df <- object@used_genesets
        
        genes <-
          gene_set_df %>%
          dplyr::filter(ont == current$gene_set) %>%
          dplyr::filter(gene %in% base::rownames(rna_assay())) %>%
          dplyr::pull(gene)
        
        if(current$method_gs == "mean"){
          
          geneset_vls <-
            base::colMeans(rna_assay()[genes, ]) %>%
            base::as.data.frame() %>%
            magrittr::set_colnames(value = "expr_score") %>%
            tibble::rownames_to_column(var = "barcodes")
          
        } else if(current$method_gs %in% c("gsva", "ssgsea", "zscore", "plage")) {
          
          shiny::showNotification(
            ui = stringr::str_c("Calculating gene set score according to method: '", current$method_gs, "'. This might take a few moments.", sep = ""),
            type = "message")
          
          geneset_vls <-
            GSVA::gsva(expr = rna_assay()[genes,], gset.idx.list = gene_set_df, mx.diff = 1, parallel.sz = 2, method = current$method_gs, verbose = F) %>%
            base::t() %>%
            base::as.data.frame() %>%
            magrittr::set_colnames(value = "expr_score") %>%
            tibble::rownames_to_column(var = "barcodes")
          
        }
        
        base::return(geneset_vls)
        
        
      })
      
      # fdata
      fdata <- shiny::reactive({
        
        fdata <-
          getFeatureDf(object = object, of_sample = current$sample)[, c("barcodes", current$feature)]
        
        base::return(fdata)
        
      })
      
      # joined data.frame
      joined_df <- shiny::reactive({
        
        if(current$color_code == "genes"){
          
          joined_df <-
            dplyr::left_join(x = sample_coords(), y = gene_vls(), by = "barcodes")
          
        } else if(current$color_code == "gene_sets"){
          
          joined_df <-
            dplyr::left_join(x = sample_coords(), y = geneset_vls(), by = "barcodes")
          
        } else if(current$color_code == "feature"){
          
          joined_df <-
            dplyr::left_join(x = sample_coords(), y = fdata(), by = c("barcodes"))
          
        }
        
        base::return(joined_df)
        
      })
      
      # variable
      variable <- shiny::reactive({
        
        if(current$color_code %in% c("genes", "gene_sets")){
          
          variable <- "expr_score"
          
        } else if(current$color_code == "feature") {
          
          variable <- current$feature
          
        }
        
        base::return(variable)
        
      })
      
      # color variable
      color_variable <- shiny::reactive({
        
        dplyr::pull(smoothed_df(), variable())
        
      })
      
      color_min <- shiny::reactive({
        
        base::min(color_variable()) %>%
          base::round(digits = 2)
        
      })
      
      color_max <- shiny::reactive({
        
        base::max(color_variable()) %>%
          base::round(digits = 2)
        
      })
      
      color_median <- shiny::reactive({
        
        stats::median(color_variable()) %>%
          base::round(digits = 2)
        
      })
      
      # smoothed_df
      smoothed_df <- shiny::reactive({
        
        shiny::validate(
          shiny::need(joined_df(), message = "Click on 'Plot & Update' to display the plot.")
        )
        
        if(base::as.numeric(input$pt_smooth) != 0){
          
          smoothed_df <-
            hlpr_smooth_shiny(coords_df = joined_df(),
                              variable = variable(),
                              smooth_span = base::as.numeric(input$pt_smooth))
          
          if(current$color_code %in% c("genes", "gene_sets")){
            
            smoothed_df <-
              purrr::imap_dfr(.x = smoothed_df,
                              .f = hlpr_normalize_imap,
                              aspect = "",
                              subset = variable())
            
          }
          
          base::return(smoothed_df)
          
        } else {
          
          if(current$color_code %in% c("genes", "gene_sets")){
            
            smoothed_df <-
              purrr::imap_dfr(.x = joined_df(),
                              .f = hlpr_normalize_imap,
                              aspect = "",
                              subset = variable())
            
            base::return(smoothed_df)
            
          } else {
            
            smoothed_df <- joined_df()
            
            base::return(smoothed_df)
            
          }
          
        }
        
      })
      
      # geom_point_add_on
      geom_point_add_on <- shiny::reactive({
        
        #color <- dplyr::pull(.data = smoothed_df(), variable())
        
        add_on <-
          list(
            ggplot2::geom_point(data = smoothed_df() %>% arrange(.data[[variable()]]),
                                mapping = ggplot2::aes(x = x, y = y, color = .data[[variable()]]),
                                size = input$pt_size,
                                alpha = (1-input$pt_alpha))
          )
        
        base::return(add_on)
        
      })
      
      #----- Scale color add-on -----#
      
      color_add_on <- shiny::reactive({
        
        color_min <- input$scale_color_min
        color_max <- input$scale_color_max
        color_mid <- input$scale_color_mid
        
        if(base::is.numeric(color_variable())){
          
          if(current$pt_clrsp %in% validColorSpectra()[["Diverging"]]){
            
            add_on <-
              confuns::scale_color_add_on(clrsp = current$pt_clrsp,
                                          limits = c(color_min,
                                                     color_max),
                                          mid = color_mid,
                                          oob = scales::squish)
            
          } else {
            
            add_on <-
              confuns::scale_color_add_on(clrsp = current$pt_clrsp,
                                          limits = c(color_min,
                                                     color_max),
                                          oob = scales::squish)
            
          }
          
        } else if(!base::is.numeric(color_variable())){
          
          add_on <-
            list(confuns::scale_color_add_on(variable = "discrete", clrp = current$pt_clrp),
                 ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))))
          
        }
        
        base::return(add_on)
        
      })
      
      
      #----- Theme add-ons -----#
      coords_add_on <- shiny::reactive({
        
        if("coords" %in% input$display_add_ons){
          
          add_on <-
            list(ggplot2::theme_bw(),
                 ggplot2::theme(
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank()
                 ))
          
        } else {
          
          add_on <-
            list(ggplot2::theme_void())
          
        }
        
        base::return(add_on)
        
      })
      
      legend_add_on <- shiny::reactive({
        
        if("legend" %in% input$display_add_ons){
          
          if(current$color_code %in% c("gene_sets", "genes")){
            
            legend_title = "Expr.\nscore"
            
          } else {
            
            legend_title = current$feature
            
          }
          
          add_on <-
            list(ggplot2::labs(color = legend_title))
          
        } else {
          
          add_on <-
            list(ggplot2::theme(legend.position = "none"))
          
        }
        
        base::return(add_on)
        
        
      })
      
      title_add_on <- shiny::reactive({
        
        if("title" %in% input$display_add_ons){
          
          if(current$color_code == "genes"){
            
            genes <- current$genes
            
            if(length(genes) > 5){
              
              genes <- c(genes[1:5], stringr::str_c("... +", (length(genes)-5), sep = " "))
              
            }
            
            genes_string <- stringr::str_c(genes, collapse = ", ")
            
            plot_title <- stringr::str_c("Genes:", genes_string, sep = " ")
            
          } else if(current$color_code == "gene_sets"){
            
            gene_set <- current$gene_set
            
            gene_set_string <- stringr::str_c(gene_set, " (", current$method_gs, ")", sep = "")
            
            plot_title <- stringr::str_c("Gene set:", gene_set_string, sep = " ")
            
          } else {
            
            plot_title <- stringr::str_c("Feature:", current$feature, sep = " ")
            
          }
          
          add_on <- ggplot2::labs(title = plot_title)
          
        } else {
          
          add_on <- NULL
          
        }
        
        base::return(add_on)
        
        
      })
      
      segmentation_add_on <- reactive({
        
        if("segmentation" %in% input$display_add_ons){
          
          if(nrow(segmentation_df()) == 0){
            
            shiny::showNotification(ui = stringr::str_c("Sample", current$sample, "has not been segmented so far.", sep = " "))
            base::return(list())
            
          } else {
            
            segm_layer <-
              list(
                ggalt::geom_encircle(data = segmentation_df(), alpha = 0.75, expand = 0.025,
                                     mapping = ggplot2::aes(x = x, y = y, group = segmentation, fill = segmentation)),
                confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = "milo", guide = FALSE)
                
              )
            
            base::return(segm_layer)
            
          }
          
        } else {
          
          base::return(list())
          
        }
        
      })
      
      segmentation_df <- reactive({
        
        segm_df <- joinWith(object = reactive_object(),
                            spata_df = getCoordsDf(reactive_object(), current$sample),
                            features = "segmentation",
                            verbose = FALSE) %>%
          dplyr::filter(!segmentation %in% c("none", ""))
        
        base::return(segm_df)
        
      })
      
      # -----
      
      # Assembled plot ----------------------------------------------------------
      
      assembled_plot <- shiny::reactive({
        
        shiny::req(input$update_plot)
        
        ggplot2::ggplot() +
          image_add_on() +
          geom_point_add_on() +
          color_add_on() +
          title_add_on() +
          segmentation_add_on() +
          ggplot2::coord_equal() +
          coords_add_on() +
          legend_add_on()
        
      })
      
      # -----
      
      # Observe events ----------------------------------------------------------
      
      # update plot by updating reactive values
      oe <- shiny::observeEvent(input$update_plot, {
        
        current$sample = input$sample_opts
        current$color_code = input$aes_clr_opts
        
        if(current$color_code == "genes"){
          
          current$genes = input$aes_clr_opts_detailed
          
        } else if(current$color_code == "gene_sets"){
          
          current$gene_set = input$aes_clr_opts_detailed
          current$method_gs = input$method_gs
          
        } else if(current$color_code == "feature"){
          
          current$feature = input$aes_clr_opts_detailed
          
        }
        
        current$pt_size = input$pt_size
        current$pt_clrsp = input$pt_clrsp
        current$pt_clrp = input$pt_clrp
        current$pt_alpha = input$pt_alpha
        
        if(base::isTRUE(input$reset_select_genes) &&
           current$color_code == "genes"){
          reset_select_genes((reset_select_genes() + 1))
        }
        
      })
      
      # -----
      
      # Return values -----------------------------------------------------------
      
      return_list <- shiny::reactive({
        
        list(
          assembled_plot = shiny::reactive({assembled_plot()}),
          dblclick = shiny::reactive({input$surface_plot_dblclick}),
          current_setting = shiny::reactive({current}),
          smoothed_df = shiny::reactive({smoothed_df()}),
          variable = shiny::reactive({variable()}),
          variable_name = shiny::reactive(input$aes_clr_opts_detailed),
          pt_size_reactive = shiny::reactive(input$pt_size)
        )
        
      })
      
      
      base::return(return_list)
      
      # -----
      
    })
  
}


#' @title Normalize gene or gene set values
#'
#' @description Helper function to use within \code{purrr::imap()}
#'
#' @param variable The variable to normalize (if matches requirements).
#' @param var_name The name of the variable to smooth.
#' @param verbose Logical
#' @param aspect Gene or Gene set
#' @param subset A character vector of variable names that are to be normalized
#' @param pb An R6 progress bar object.
#'
#' @return A normalized variable (data.frame within \code{purrr::imap()})
#' @export

hlpr_normalize_imap <- function(variable,
                                var_name,
                                aspect,
                                subset){
  
  if(!base::is.numeric(variable) | !var_name %in% subset){
    
    base::return(variable)
    
  } else if(base::all(variable == 0)){
    
    if(var_name == "mean_genes"){
      
      var_name <- "average"
      
    }
    
    base::warning(stringr::str_c(aspect, var_name, "contains only 0s. Returning NULL.", sep = " "))
    base::return(NULL)
    
  } else if(base::length(base::unique(variable)) == 1){
    
    if(var_name == "mean_genes"){
      
      var_name <- "average"
      
    }
    
    base::warning(stringr::str_c(aspect, var_name, "is uniformly expressed. Returning NULL.", sep = " "))
    base::return(NULL)
    
  } else {
    
    # normalize variable
    res <-
      (variable - base::min(variable)) /
      (base::max(variable) - base::min(variable))
    
    if(!base::any(base::is.na(res))){
      
      base::return(res)
      
    } else {
      
      base::warning(stringr::str_c(aspect, var_name, "normalization resulted in NaNs. Returning NULL.", sep = " "))
      base::return(NULL)
      
    }
    
  }
  
}


#' @title Smooth variable spatially in mini-shiny-apps
#'
#' @description Helper function to use independently (or in a pipe)
#'
#' @inherit hlpr_smooth params
#'
#' @return Data.frame with the smoothed variable specified in \code{variable}.
#'
#'@export
hlpr_smooth_shiny <- function(variable,
                              coords_df,
                              smooth_span){
  
  base::colnames(coords_df)[base::which(base::colnames(coords_df) == variable)] <- "response_variable"
  
  if(base::is.numeric(coords_df$response_variable)){
    
    model <- stats::loess(formula = response_variable ~ x * y, span = smooth_span, data = coords_df)
    
    smoothed_df_prel <-
      broom::augment(model) %>%
      dplyr::select(x, y, .fitted) %>%
      magrittr::set_colnames(value = c("x", "y", variable))
    
    selected_df <- dplyr::select(coords_df, -c("x", "y", "response_variable"))
    
    smoothed_df <-
      base::cbind(smoothed_df_prel, selected_df) %>%
      dplyr::select(barcodes, sample, x, y, dplyr::everything()) %>%
      as.data.frame()
    
    
    # if coords_df derived from trajectory analysis
    if("trajectory_order" %in% base::colnames(coords_df)){
      
      smoothed_df$trajectory_order <- coords_df$trajectory_order
      
    }
    
    if(base::nrow(smoothed_df) == base::nrow(coords_df)){
      
      return(smoothed_df)
      
    } else {
      
      shiny::showNotification(ui = "Smoothing failed. Return original values.",
                              type = "warning")
      
      return(coords_df)
      
    }
    
    
  } else {
    
    shiny::showNotification(ui = "Can not smooth features that aren't of class 'numeric'. Skip smoothing.",
                            type = "warning")
    
    base::colnames(coords_df)[base::which(base::colnames(coords_df) == "response_variable")] <- variable
    
    return(coords_df)
    
  }
  
}



