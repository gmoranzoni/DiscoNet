#' Visualize Network
#'
#' Creates static and interactive visualizations of a given network with customizable layouts and aesthetic attributes.
#' It allows the user to visualize networks using discrete or continuous data for node color and provides several layout options.
#'
#' @param igraph The network in the form of an `igraph` object.
#' @param layout The desired network layout, with default being `with_fr()`. Supported layouts include `with_fr()`, `with_dh()`, `with_gem()`,
#' `with_graphopt()`, `with_kk()`, `with_lgl()`, `with_mds()`, `as_star()`, `as_tree()`, `in_circle()`, `nicely()`, `on_grid()`, and `randomly()`.
#' @param discrete The name of the column storing the discrete values that should be used for node color mapping. Default is `NULL`.
#' @param continuous The name of the column storing the continuous variable that should be used for node color mapping. Default is `NULL`.
#'
#' @return A list containing the static ggplot object and interactive plotly visualization of the network.
#'
#' @import igraph
#' @import ggplot2
#' @importFrom ggnetwork ggnetwork
#' @importFrom dplyr as_tibble %>%
#' @importFrom ggplot2 ggplot scale_fill_viridis_d scale_fill_discrete scale_color_viridis_c scale_shape_manual scale_color_continuous labs
#' @importFrom ggnetwork geom_edges geom_nodes
#' @importFrom ggraph theme_graph
#' @importFrom plotly ggplotly style
#' @importFrom forcats as_factor
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming `sample_igraph` is a predefined igraph object, `sample_layout` a chosen layout,
#' # and `sample_discrete` & `sample_continuous` are column names in the igraph object.
#' visualize_network(igraph = sample_igraph, layout = sample_layout, discrete = sample_discrete, continuous = sample_continuous)
#'}
visualize_network <- function(igraph, layout = with_fr(), discrete = NULL, continuous = NULL) {

  # Store the original column names before they get overwritten
  discrete_colname <- discrete
  continuous_colname <- continuous

  # remove duplicate edges
  igraph <- simplify(igraph, remove.multiple = TRUE)

  # preparation of the network
  df_net <- ggnetwork(igraph, layout = layout)
  df_net <- df_net %>% as_tibble()

  # conversion of the column seed from numeric to factor.
  df_net$seed <- as_factor(as.numeric(df_net$seed))

  # extract the discrete and continuous variable from df_net, to use them in the plot
  if(!is.null(discrete)) {
    discrete <- df_net[[discrete]]
  }
  if(!is.null(continuous)) {
    continuous <- df_net[[continuous]]
  }

  if(!is.null(discrete) & !is.null(continuous)) {

    plot <- df_net %>%
      ggplot(., aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(linewidth = 0.4, alpha = 0.25) +
      geom_nodes(aes(text = name, fill = .data[[discrete_colname]], color = .data[[continuous_colname]], shape = seed), size = 4, stroke = 1) +
      scale_fill_viridis_d(option = "plasma", begin = 0.3, end = 0.9, na.value = "#6189B8") +
      # scale_fill_discrete(na.value = "#6189B8", type = c("#F4157D","#FFF70A")) +
      scale_color_viridis_c(option = "viridis", na.value = "#6189B8") + # apply viridis colors
      #scale_color_continuous(low = "#FF4365", high = "#9FD356", na.value = "#6189B8") +
      scale_shape_manual(values = c("0" = 21, "1" = 24)) +
      theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black') +
      labs()


  } else if(!is.null(discrete) & is.null(continuous)) {
    plot <- df_net %>%
      ggplot(., aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(linewidth = 0.4, alpha = 0.25) +
      geom_nodes(aes(text = name, fill = .data[[discrete_colname]], shape = seed), color = "blue", size = 4, stroke = 1) +
      scale_fill_discrete(na.value = "#6189B8", type = c("#F4157D","#FFF70A")) +
      scale_color_viridis(option = "viridis", na.value = "#6189B8") + # apply viridis colors
      #scale_color_continuous(guide = "none") +
      scale_shape_manual(values = c("0" = 21, "1" = 24)) +
      theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black') +
      labs()
  } else if(is.null(discrete) & !is.null(continuous)) {
    plot <- df_net %>%
      ggplot(., aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(linewidth = 0.4, alpha = 0.25) +
      geom_nodes(aes(text = name, color = .data[[continuous_colname]], shape = seed), fill = "red", size = 4, stroke = 1) +
      scale_fill_discrete(guide = "none") +
      scale_color_viridis(option = "viridis", na.value = "#6189B8") + # apply viridis colors
      #scale_color_continuous(low = "#FF4365", high ="#9FD356", na.value = "#6189B8") +
      scale_shape_manual(values = c("0" = 21, "1" = 24)) +
      theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black') +
      labs()

  } else if(is.null(discrete) & is.null(continuous)) {
    plot <- df_net %>%
      ggplot(., aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(linewidth = 0.4, alpha = 0.25) +
      geom_nodes(aes(text = name, shape = seed), color = "#6189B8", fill = "#6189B8", size = 4, stroke = 1) +
      scale_fill_discrete(guide = "none") +
      scale_color_continuous(guide = "none") +
      scale_shape_manual(values = c("0" = 21, "1" = 24)) +
      theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black') +
      labs()
  }

  # Get the column names for tooltips
  tooltip_names <- c("text", "seed")
  if (!is.null(discrete_colname)) {
    tooltip_names <- c(tooltip_names, discrete_colname)
  }
  if (!is.null(continuous_colname)) {
    tooltip_names <- c(tooltip_names, continuous_colname)
  }

  font <- list(
    color = "#6189B8",
    size = 16,
    family = "AppleGothic")

  label <- list(
    bgcolor = "#CAD5E3",
    bordercolor = "#6189B8",
    font = font,
    hoverformat = ".2f")


  plot_plotly <- plot %>%
    ggplotly(tooltip = tooltip_names) %>%
    style(hoverlabel = label)

  plot_list <- list()
  plot_list$static <- plot
  plot_list$interactive <- plot_plotly

  return(plot_list)
}
