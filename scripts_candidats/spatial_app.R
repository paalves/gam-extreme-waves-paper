
rm(list=ls())

library(shiny)
library(leaflet)
library(leaflet.extras) 
library(leaflet.minicharts)
library(resourcecode)
library(resourcecodedata)
library(dplyr)
library(ggplot2)
library(sf)
library(interp) 
library(tibble)
library(lubridate)
library(tidyr)
library(RColorBrewer)

# ==============================================================================
# 0. CONFIGURATION & DEFINITIONS
# ==============================================================================

VARS_AVAILABLE_FLAT <- c(
  "Hauteur significative (m)" = "hs",
  "Période du pic (s)" = "tp",
  "Cambrure" = "steepness",
  "Direction vagues (°)" = "dp",
  
  "Vitesse courant (m/s)" = "Cspd",
  "Direction courant (°)" = "Cdirdeg",
  "Hauteur d'eau (m)" = "wlv",
  "Ecart angulaire vague-courant (°)" = "angleDiff",
  "Profil vagues-courant (catégorie)" = "wcProfile", 
  
  "Vitesse vent (m/s)" = "Wspd",
  "Direction vent (°)" = "Wdirdeg",
  
  "Cycle annuel (0-1)" = "cycle"
)

# Variables considered "Directional" for Polar Plots
DIRECTIONAL_VARS <- c("dp", "Cdirdeg", "Wdirdeg", "angleDiff")

# Function to get Label from Var Name
get_label <- function(var_name) {
  lbl <- names(VARS_AVAILABLE_FLAT)[VARS_AVAILABLE_FLAT == var_name]
  if(length(lbl) == 0) return(var_name)
  return(lbl)
}

transform_dataset <- function(data) {
  df <- data %>%
    mutate(mois=month(time), annee=year(time), jour=day(time), heure=hour(time))
  
  if(all(c("uubr", "vubr") %in% names(df))) df <- df %>% mutate(OrbSpeed = sqrt(uubr^2 + vubr^2))
  if(all(c("hs", "tp") %in% names(df))) df <- df %>% mutate(steepness = hs / tp)
  
  df <- df %>%
    mutate(
      current_datetime = as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00", annee, mois, jour, heure), tz = "UTC"),
      year_start = as.POSIXct(sprintf("%04d-01-01 00:00:00", annee), tz = "UTC"),
      cycle = as.numeric(difftime(current_datetime, year_start, units = "hours"))
    ) %>%
    group_by(annee) %>%
    mutate(cycle = cycle / max(cycle, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(-current_datetime, -year_start)
  
  if(all(c("ucur", "vcur") %in% names(df))) {
    df <- df %>%
      mutate(
        Cspd = sqrt(ucur^2 + vcur^2), 
        Cdir = atan2(vcur, ucur), 
        Cdirdeg = (90 - Cdir * (180/pi)) %% 360
      )
  }
  
  if(all(c("uwnd", "vwnd") %in% names(df))) df <- df %>% mutate(Wspd = sqrt(uwnd^2 + vwnd^2), Wdir = atan2(vwnd, uwnd), Wdirdeg = (90 - Wdir * (180/pi)) %% 360)
  
  if(all(c("dp", "Cdirdeg") %in% names(df))) {
    df <- df %>%
      mutate(
        dpAngleDiff = case_when(dp>0 & dp<180 ~ dp+180, dp>180 & dp<360 ~ dp-180, dp==0 | dp==360 ~ 180, TRUE ~ 0),
        angleDiff = (dpAngleDiff - Cdirdeg) %% 360
      ) %>%
      mutate(
        wcProfile = case_when(
          angleDiff >= 150 & angleDiff <= 210 ~ "Contre-courant",
          angleDiff >= 330 | angleDiff <= 30  ~ "Co-courant",
          TRUE ~ "Perpendiculaire"
        )
      )
  }
  
  return(df)
}

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(
  navbarPage("Statistiques spatiales Resourcecode",
    
    # --- ONGLET 1 : IMPORTATION ---
    tabPanel("Importation des données",
      sidebarLayout(
        sidebarPanel(
          h3("1. Définir la zone"),
          helpText("Dessinez un rectangle sur la carte."),
          verbatimTextOutput("zone_info"),
          
          hr(),
          
          h3("2. Configuration"),
          numericInput("limit_points", "Nombre de points à importer :", value=50, min=10, max=500),
          dateRangeInput("date_range", "Période d'analyse :", 
                         start="1994-01-01", end="1994-01-08", 
                         min="1994-01-01", max="2020-12-31"),
          
          hr(),
          h4("Variables à traiter"),
          checkboxGroupInput("selected_families", "Familles de variables :",
                             choices = c("Vagues", "Courants", "Vent"),
                             selected = c("Vagues", "Courants", "Vent")),
          
          hr(),
          actionButton("run_analysis", "Importation des données", class="btn-primary", width="100%"),
          br(), br(),
          textOutput("status_text")
        ),
        mainPanel(
          leafletOutput("map_select", height = "700px")
        )
      )
    ),
    
    # --- ONGLET 2 : ANALYSE SPATIALE ---
    tabPanel("Analyse spatiale",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Visualisation spatiale"),
          selectInput("var_to_plot", "Variable à afficher :", choices = NULL),
          radioButtons("stat_type", "Statistique :", 
                       choices=c("Moyenne"="mean", "Maximum"="max", "Ecart-type"="sd", "Minimum"="min")),
          hr(),
          checkboxInput("use_interp", "Interpolation", value = TRUE),
          conditionalPanel("input.use_interp == false",
             sliderInput("point_size", "Taille des points :", min=0.5, max=10, value=3, step=0.5)
          )
        ),
        mainPanel(
          width = 9,
          plotOutput("spatial_plot", height = "700px")
        )
      )
    ),
    
    # --- ONGLET 3 : MINIGRAPHIQUES ---
    tabPanel("Mini-graphiques",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Minicharts"),
          
          selectInput("mc_var", "Variable à afficher (Y) :", choices = NULL),
          selectInput("mc_filter_var", "Variable de filtrage (X) :", choices = NULL),
          selectInput("mc_type", "Type de graph :", 
                      choices = c("Histogramme / barplot" = "bar")),
          
          conditionalPanel("input.mc_type == 'polar'", helpText("Rose")),
          conditionalPanel("input.mc_type == 'bar'",
             helpText("1. On garde les données où la variable X est > quantile."),
             helpText("2. On affiche l'histogramme de la variable Y."),
             sliderInput("mc_quantile", "Filtre quantile local (par noeud) :", min=0, max=0.99, value=0, step=0.05),
             sliderInput("mc_bins", "Nombre de barres :", min=3, max=10, value=5)
          ),
          
          hr(),
          actionButton("update_mc", "Mettre à jour la carte", class="btn-success")
        ),
        mainPanel(
          width = 9,
          leafletOutput("map_minicharts", height = "700px")
        )
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    bbox = NULL,           
    selected_nodes = NULL,
    data_res = NULL,
    data_ts = NULL,
  )
  
  # --------------------------------------------------------------------------
  # 1. Selection & Zone (Onglet 1)
  # --------------------------------------------------------------------------
  observeEvent(input$map_select_draw_new_feature, {
    feat <- input$map_select_draw_new_feature
    coords <- feat$geometry$coordinates[[1]]
    lons <- sapply(coords, function(x) x[[1]])
    lats <- sapply(coords, function(x) x[[2]])
    rv$bbox <- list(xmin = min(lons), xmax = max(lons), ymin = min(lats), ymax = max(lats))
    
    in_zone <- rscd_field %>%
      filter(longitude >= rv$bbox$xmin, longitude <= rv$bbox$xmax,
             latitude >= rv$bbox$ymin, latitude <= rv$bbox$ymax)
    rv$selected_nodes <- in_zone
  })
  
  output$zone_info <- renderText({
    if(is.null(rv$selected_nodes)) return("Aucune zone sélectionnée.")
    paste("Points sélectionnés :", nrow(rv$selected_nodes))
  })
  
  # --------------------------------------------------------------------------
  # 2. Global Calculation (Onglet 1)
  # --------------------------------------------------------------------------
  observeEvent(input$run_analysis, {
    req(rv$selected_nodes)
    nodes_df <- rv$selected_nodes
    if(nrow(nodes_df) == 0) { showNotification("Zone vide.", type="error"); return() }
    
    if(nrow(nodes_df) > input$limit_points) {
      nodes_to_process <- nodes_df %>% sample_n(input$limit_points)
      showNotification(paste("Echantillonnage :", input$limit_points, "points."), type="warning")
    } else {
      nodes_to_process <- nodes_df
    }
    
    params_to_fetch <- c()
    if("Vagues" %in% input$selected_families) params_to_fetch <- c(params_to_fetch, "hs", "tp", "dp")
    if("Courants" %in% input$selected_families) params_to_fetch <- c(params_to_fetch, "ucur", "vcur", "wlv")
    if("Vent" %in% input$selected_families) params_to_fetch <- c(params_to_fetch, "uwnd", "vwnd")
    
    if(length(params_to_fetch) == 0) {
      showNotification("Veuillez sélectionner au moins une famille.", type="error")
      return()
    }
    
    output$status_text <- renderText("Récupération des données...")
    progress <- shiny::Progress$new()
    progress$set(message = "Importation en cours :", value = 0)
    on.exit(progress$close())
    
    n_iter <- nrow(nodes_to_process)
    res_stats_list <- list()
    res_ts_list <- list()
    
    for(i in 1:n_iter) {
      if(i %% 5 == 0) progress$set(value = i/n_iter, detail = paste(i, "/", n_iter))
      
      tryCatch({
        df_raw <- get_parameters(
          node = nodes_to_process$node[i],
          start = as.POSIXct(paste0(input$date_range[1], " 00:00:00"), tz="UTC"),
          end   = as.POSIXct(paste0(input$date_range[2], " 23:00:00"), tz="UTC"),
          parameters = params_to_fetch
        )
        
        df_transformed <- transform_dataset(df_raw)
        
        df_transformed$node <- nodes_to_process$node[i]
        df_transformed$longitude <- nodes_to_process$longitude[i]
        df_transformed$latitude <- nodes_to_process$latitude[i]
        
        res_ts_list[[i]] <- df_transformed
        
        stats_row <- df_transformed %>%
          select(where(is.numeric)) %>%
          select(-any_of(c("node", "latitude", "longitude", "annee", "mois", "jour", "heure", "time"))) %>%
          summarise(across(everything(), list(
            mean = ~mean(.x, na.rm=TRUE),
            max = ~max(.x, na.rm=TRUE),
            min = ~min(.x, na.rm=TRUE),
            sd = ~sd(.x, na.rm=TRUE)
          ), .names = "{.col}_{.fn}"))
        
        stats_row$node <- nodes_to_process$node[i]
        stats_row$longitude <- nodes_to_process$longitude[i]
        stats_row$latitude <- nodes_to_process$latitude[i]
        
        res_stats_list[[i]] <- stats_row
        
      }, error = function(e) {})
    }
    
    if(length(res_stats_list) > 0) {
      rv$data_res <- bind_rows(res_stats_list)
      rv$data_ts <- bind_rows(res_ts_list)
      
      cols_present <- names(rv$data_ts)
      
      valid_vars_display <- VARS_AVAILABLE_FLAT[VARS_AVAILABLE_FLAT %in% cols_present]
      numeric_cols <- rv$data_ts %>% select(where(is.numeric)) %>% names()
      valid_vars_filter <- valid_vars_display[valid_vars_display %in% numeric_cols]
      
      updateSelectInput(session, "var_to_plot", choices = valid_vars_display)
      updateSelectInput(session, "mc_var", choices = valid_vars_display)
      updateSelectInput(session, "mc_filter_var", choices = valid_vars_filter, selected = valid_vars_filter[1])
      
      output$status_text <- renderText(paste("Terminé ! Points importés :", nrow(rv$data_res)))
    } else {
      output$status_text <- renderText("Erreur: aucune donnée récupérée.")
    }
  })
  
  # --------------------------------------------------------------------------
  # 3. Dynamic UI Updates (Onglet 3)
  # --------------------------------------------------------------------------
  observeEvent(input$mc_var, {
    req(input$mc_var)
    var_code <- input$mc_var
    is_directional <- var_code %in% DIRECTIONAL_VARS
    
    if(is_directional) {
      choices_mc <- c("Histogramme / barplot" = "bar", "Rose" = "polar")
    } else {
      choices_mc <- c("Histogramme / barplot" = "bar")
    }
    updateSelectInput(session, "mc_type", choices = choices_mc)
  })
  
  # --------------------------------------------------------------------------
  # 4. MAPS (Onglet 1 & 2)
  # --------------------------------------------------------------------------
  output$map_select <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE)) %>% 
      addTiles() %>%
      setView(lng = -4, lat = 48, zoom = 6) %>%
      addDrawToolbar(
        targetGroup = "draw",
        rectangleOptions = drawRectangleOptions(shapeOptions = drawShapeOptions(fillOpacity = 0.2, color="red")),
        polylineOptions = FALSE, polygonOptions = FALSE, 
        circleOptions = FALSE, markerOptions = FALSE, circleMarkerOptions = FALSE,
        editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions())
      ) %>%
      addCircleMarkers(data = rscd_field, 
                       lng=~longitude, lat=~latitude, 
                       radius=2, color="black", stroke=FALSE, fillOpacity=0.5, 
                       group="percu")
  })
  
  output$spatial_plot <- renderPlot({
    req(rv$data_res, rv$bbox, input$var_to_plot)
    
    df_plot_full <- rv$data_res
    target_col <- paste0(input$var_to_plot, "_", input$stat_type)
    
    if(!target_col %in% names(df_plot_full)) return(ggplot() + annotate("text", x=0, y=0, label="Donnée/stat non dispo"))
    
    df_viz <- df_plot_full %>% select(longitude, latitude, z = all_of(target_col))
    
    x_margin <- (rv$bbox$xmax - rv$bbox$xmin) * 0.1
    y_margin <- (rv$bbox$ymax - rv$bbox$ymin) * 0.1
    x_lims <- c(rv$bbox$xmin - x_margin, rv$bbox$xmax + x_margin)
    y_lims <- c(rv$bbox$ymin - y_margin, rv$bbox$ymax + y_margin)
    
    p <- ggplot()
    
    if(input$use_interp) {
      in_view_mask <- rscd_field$longitude >= x_lims[1] & rscd_field$longitude <= x_lims[2] &
                      rscd_field$latitude >= y_lims[1] & rscd_field$latitude <= y_lims[2]
      
      t1 <- in_view_mask[resourcecodedata::rscd_triangles[1,]]
      t2 <- in_view_mask[resourcecodedata::rscd_triangles[2,]]
      t3 <- in_view_mask[resourcecodedata::rscd_triangles[3,]]
      tri_keep <- t1 & t2 & t3
      subset_tri <- resourcecodedata::rscd_triangles[, tri_keep, drop=FALSE]
      
      if(ncol(subset_tri) > 0 && nrow(df_viz) > 5) {
        needed_node_ids <- unique(as.vector(subset_tri))
        target_nodes <- rscd_field[needed_node_ids, ]
        
        tryCatch({
          interp_res <- interp::interp(
            x = df_viz$longitude, y = df_viz$latitude, z = df_viz$z,
            xo = target_nodes$longitude, yo = target_nodes$latitude,
            output = "points", linear = FALSE, extrap = TRUE
          )
          z_full <- rep(NA, nrow(rscd_field))
          z_full[needed_node_ids] <- interp_res$z
          
          xyzgz <- tibble::tibble(
            x = rscd_field$longitude[as.vector(subset_tri)],
            y = rscd_field$latitude[as.vector(subset_tri)],
            z = z_full[as.vector(subset_tri)],
            g = rep(1:ncol(subset_tri), each=3)
          ) %>% group_by(g) %>% filter(!any(is.na(z))) %>% ungroup()
          
          p <- p + geom_polygon(data = xyzgz, aes(x, y, group = g, fill = z, col = z)) +
            scale_fill_viridis_c(option="C", name=get_label(input$var_to_plot)) +
            scale_color_viridis_c(option="C", guide="none")
        }, error = function(e) {
          p <<- p + geom_point(data=df_viz, aes(x=longitude, y=latitude, color=z), size=3)
        })
      } else {
        p <- p + geom_point(data=df_viz, aes(x=longitude, y=latitude, color=z), size=3)
      }
    } else {
      p <- p + geom_point(data=df_viz, aes(x=longitude, y=latitude, color=z), 
                          size=input$point_size, shape=16) +
        scale_color_viridis_c(option="C", name=get_label(input$var_to_plot))
    }
    
    coast_vis <- rscd_coastline %>% 
      filter(longitude >= x_lims[1], longitude <= x_lims[2], latitude >= y_lims[1], latitude <= y_lims[2])
    if(nrow(coast_vis) > 0) p <- p + geom_path(data=coast_vis, aes(x=longitude, y=latitude), color="black", linewidth=0.3)
    
    islands_vis <- rscd_islands %>%
      filter(longitude >= x_lims[1], longitude <= x_lims[2], latitude >= y_lims[1], latitude <= y_lims[2])
    if(nrow(islands_vis) > 0) p <- p + geom_path(data=islands_vis, aes(x=longitude, y=latitude, group=ID), color="black", linewidth=0.3)
    
    p + coord_sf(xlim=x_lims, ylim=y_lims, expand=FALSE, crs=sf::st_crs(4326)) +
      theme_void() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "right") +
      labs(title = get_label(input$var_to_plot), subtitle = paste("Statistique :", input$stat_type))
  })
  
  output$map_minicharts <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      setView(lng = -4, lat = 48, zoom = 6)
  })
  
  # --------------------------------------------------------------------------
  # 5. MINICHARTS LOGIC (Onglet 3)
  # --------------------------------------------------------------------------
  observeEvent(input$update_mc, {
    req(rv$data_ts, input$mc_var, input$mc_filter_var)
    
    df_full <- rv$data_ts
    
    var_display <- input$mc_var
    var_filter  <- input$mc_filter_var
    
    df_work <- df_full %>%
      select(node, longitude, latitude, time, 
             val = all_of(var_display), 
             filter_val = all_of(var_filter)) %>%
      filter(!is.na(val), !is.na(filter_val))
    
    if(nrow(df_work) == 0) return()
    
    df_filt <- df_work %>%
      group_by(node) %>%
      mutate(threshold = quantile(filter_val, input$mc_quantile, na.rm = TRUE)) %>%
      filter(filter_val >= threshold) %>%
      ungroup()
    
    if(nrow(df_filt) == 0) {
      showNotification("Données vides après filtrage.", type="warning")
      return()
    }
    
    is_numeric <- is.numeric(df_filt$val)
    mat_data <- NULL
    colors_vec <- NULL
    chart_type <- "bar"
    w_size <- 100 
    h_size <- 100
    
    if(input$mc_type == "polar" && is_numeric) {
      
      chart_type <- "polar-area"
      w_size <- 120 
      h_size <- 120
      
      df_filt <- df_filt %>% mutate(val = val %% 360) 
      df_filt$sector <- cut(
        (df_filt$val + 22.5) %% 360,
        breaks = seq(0, 360, 45),
        labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
        include.lowest = TRUE
      )
      
      chart_data_wide <- df_filt %>%
        count(node, longitude, latitude, sector) %>%
        pivot_wider(names_from = sector, values_from = n, values_fill = 0) %>%
        arrange(node)
      
      all_dirs <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
      for(d in all_dirs) {
        if(!d %in% names(chart_data_wide)) chart_data_wide[[d]] <- 0
      }
      chart_data_wide <- chart_data_wide %>% select(node, longitude, latitude, all_of(all_dirs))
      mat_data <- chart_data_wide %>% select(-node, -longitude, -latitude) %>% as.matrix()
      colors_vec <- brewer.pal(8, "Set2")
      
    } else {
      
      if(is_numeric) {
        
        breaks <- pretty(df_filt$val, n = input$mc_bins)
        df_filt$bin <- cut(df_filt$val, breaks = breaks, include.lowest = TRUE)
        
        chart_data_wide <- df_filt %>%
          count(node, longitude, latitude, bin) %>%
          pivot_wider(names_from = bin, values_from = n, values_fill = 0) %>%
          arrange(node)
        
        mat_data <- chart_data_wide %>% select(-node, -longitude, -latitude) %>% as.matrix()
        colors_vec <- colorRampPalette(c("#fbd9d3", "#900000"))(ncol(mat_data))
        
      } else {
        if(var_display == "wcProfile") {
          needed_cols <- c("Co-courant", "Contre-courant", "Perpendiculaire")
          chart_data_wide <- df_filt %>%
            count(node, longitude, latitude, val) %>%
            pivot_wider(names_from = val, values_from = n, values_fill = 0)
          
          for(col in needed_cols) { if(!col %in% names(chart_data_wide)) chart_data_wide[[col]] <- 0 }
          
          chart_data_wide <- chart_data_wide %>% select(node, longitude, latitude, all_of(needed_cols)) %>% arrange(node)
          mat_data <- as.matrix(chart_data_wide[, needed_cols])
          colors_vec <- c("green", "red", "yellow") 
        } else {
          chart_data_wide <- df_filt %>%
            count(node, longitude, latitude, val) %>%
            pivot_wider(names_from = val, values_from = n, values_fill = 0) %>%
            arrange(node)
          mat_data <- chart_data_wide %>% select(-node, -longitude, -latitude) %>% as.matrix()
          colors_vec <- colorRampPalette(brewer.pal(min(8, max(3, ncol(mat_data))), "Set1"))(ncol(mat_data))
        }
      }
    }
    
    leafletProxy("map_minicharts") %>%
      clearGroup("mc_anchors") %>% 
      addCircleMarkers(
        lng = chart_data_wide$longitude,
        lat = chart_data_wide$latitude,
        radius = 5,
        color = "red", 
        fillColor = "red",
        stroke = FALSE,
        fillOpacity = 1,
        group = "mc_anchors",
        popup = paste("Node:", chart_data_wide$node)
      ) %>%
      clearMinicharts() %>%
      addMinicharts(
        lng = chart_data_wide$longitude,
        lat = chart_data_wide$latitude,
        type = chart_type,
        chartdata = mat_data,
        colorPalette = colors_vec,
        width = w_size, height = h_size,
        legend = TRUE
      )
  })
}

shinyApp(ui, server)
