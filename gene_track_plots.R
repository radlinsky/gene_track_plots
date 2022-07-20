#' Plot UCSC-like tracks
#'
#' @param track_tb tibble of features minimally including:
#' @param start: exon/feature region start
#' @param end: exon/feature region end
#' @param seqid: What to label X-axis (e.g., "Chromosome 1")
#' @param add_text name of column to label feature regions, e.g., exon number (default: NA)
#' @param facet_wrap_by name of column to facet wrap the features, e.g., transcript id (default: NA)
#' @param strip_position position of the facet strip labels (default: "right")
#' @param strip_text_angle angle of the strip labels (default: 0)
#' @param xlim xlimits are automatically configured unless length 2 vector provided here (default: NA)
#' @param track_breathing_room_factor white space above/below features (default: 1)
#' 1 means track takes up all the space per facet
#' 2 means track takes up 1/2 of the space per facet..
#' 3 means track takes up 1/3 of the space per facet..
#' @param strand which strand is your main feature of interest?  (default: "+")
#' if "-" strand, the genomic coordinates will appropriately decrease left to right.
#'  Useful if negative strand gene to orient it 5' to 3' left to right...
#' @param exon_height relative height of exons given (default: 0.25)
#' @param add_connections use `jstart` and `jend` columns to add introns between exons?
#' @param color_connections Boolean. If TRUE, color the junctions as an aes
#' by creating column 'junction' = paste0(jstart, '-', jend). Default=FALSE.
#' @param intron_width width of intron lines
#' @param feature_outline_width width of exon box outlines
#'
#' @return
#' @export
#'
#' @examples
plot_tracks <- function(
    track_tb,
    strand,
    seqid,
    add_text=NA,
    facet_wrap_by=NA,
    strip_position="right",
    strip_text_angle=0,
    xlim=NA,
    track_breathing_room_factor = 1,
    exon_height=0.25,
    add_connections = FALSE,
    color_tracks = FALSE,
    base_font_size = 10,
    intron_width = 2,
    feature_outline_width = 0.25,
    exon_label_repel = FALSE,
    color_connections = FALSE
) {
    theme_size <- base_font_size
    # crude scale to match font and geom_text size
    geom_text_size <- theme_size / (14/5)
    # theme.size = (14/5) * geom.text.size
    if (add_connections) {
        track_tb <- add_introns(track_tb, groupby = facet_wrap_by)
    }
    # if introns are here, we can color them
    if (nrow(dplyr::filter(track_tb, !is.na(jstart) & !is.na(jend))) > 0) {
        # unless we dont want colored introns
        if (!color_connections) {
            plot <-
                dplyr::filter(track_tb, !is.na(jstart) & !is.na(jend)) %>%
                ggplot2::ggplot() +
                ggplot2::geom_segment(ggplot2::aes(
                    x = jstart, xend = jend,
                    y = 0, yend = 0
                ),
                size = intron_width
            )
        } else{
            plot <-
                dplyr::filter(track_tb, !is.na(jstart) & !is.na(jend)) %>%
                dplyr::mutate(junction = paste0(jstart, "-", jend)) %>%
                ggplot2::ggplot() +
                ggplot2::geom_segment(ggplot2::aes(
                    x = jstart, xend = jend,
                    y = 0, yend = 0,
                    color = junction
                ),
                size = intron_width
            )
        }

    } else{
        plot <- ggplot2::ggplot()
    }


    # if coloring tracks, fill by facet_wrap
    if (color_tracks) {
        if (length(facet_wrap_by) > 1){
            stop("Only one variable supported for aes(fill). (need to only facet_wrap_by one column)")
        }
        plot <-
            plot +
            ggplot2::geom_rect(
                data = track_tb,
                ggplot2::aes(
                    xmin = start,
                    xmax = end,
                    ymin = -exon_height / track_breathing_room_factor,
                    ymax = exon_height / track_breathing_room_factor,
                    fill = !!sym(facet_wrap_by)
                ),
                colour = "black",
                size = feature_outline_width
            )
    } else {
        # else not coloring tracks, make sure they are grey
        plot <-
            plot +
            ggplot2::geom_rect(
                data = track_tb,
                ggplot2::aes(
                    xmin = start,
                    xmax = end,
                    ymin = -exon_height / track_breathing_room_factor,
                    ymax = exon_height / track_breathing_room_factor,
                ),
                # fill = ggpubr::get_palette("Greys", 10)[4],
                fill = "#C6C6C6",
                colour = "black",
                size = feature_outline_width
            )
    # plot + scale_fill_identity()
    }

    if (exon_label_repel) {
        jitter_y <- 0.1 * exon_height
        jitter_x <- 0
    } else{
        jitter_y <- 0
        jitter_x <- 0
    }
    if (!is.na(add_text)) {
        plot <-
            plot +
            ggplot2::geom_text(
                aes(
                    x = (start + end) / 2,
                    y = # (exon_height) *
                        # (1 + 0.1), # put it on top
                        0, # put in in the middle
                    label = !!sym(add_text)
                ),
                data = track_tb,
                size = geom_text_size,
                position = position_jitter(
                    width = jitter_x,
                    height = jitter_y
                )
            )
        # nvm, was 0.25
        scale_y_by <- 0
    } else{
        scale_y_by <- 0
    }
    if (!is.na(facet_wrap_by) & strip_position != "none") {
        plot <-
            plot +
            ggplot2::facet_wrap(as.formula(paste("~", facet_wrap_by)),
                                ncol = 1,
                                strip.position = strip_position)
    }

    plot <-
        plot +
        ggplot2::scale_y_continuous(
            limits = c(
                # how much breathing room above/below track
                # scale_y_by shifts it down a bit if adding exon numbers for e.g.
                -exon_height * (1 - scale_y_by),
                exon_height * (1 + scale_y_by)
            )
        ) +
        # ggpubr::theme_pubclean(flip = TRUE, base_size = base_font_size)
        my_theme_pubclean(flip = TRUE, base_size = theme_size)

    # # keep X, Y capitalized
    # chrom <- gsub("chr", "", seqid)
    # chrom <- gsub("Chr", "", chrom)
    # chrom <- gsub("CHR", "", chrom)

    if (is.na(xlim[1]) & strand == "+") {
        plot <-
            plot +
            ggplot2::scale_x_continuous(name = seqid)
        # else flip genomic coords to orient 5' to 3' properly
    } else if (is.na(xlim[1]) & strand == "-") {
        plot <-
            plot +
            ggplot2::scale_x_reverse(name = seqid)
    } else if (strand == "+") {
        plot <-
            plot +
            ggplot2::scale_x_continuous(
                name = seqid,
                limits = xlim
            )
        # else - strand so flip genomic coords to orient 5' to 3' properly
    } else {
        plot <-
            plot +
            ggplot2::scale_x_reverse(
                name = seqid,
                limits = xlim
            )
    }

    if (strip_position == "top") {
        plot <- plot +
            ggplot2::theme(
                axis.line.y = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank()
            )
    } else if (strip_position == "right") {
        plot <- plot +
            ggplot2::theme(
                axis.line.y = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                # old (current gtx version) ggplot way to define strip angle
                strip.text.y = element_text(angle = strip_text_angle),
                # latest ggplot allows for this:
                # strip.text.y.right = element_text(angle = strip_text_angle)
            )
    } else if (strip_position == "none") {
        plot <- plot +
            ggplot2::theme(
                axis.line.y = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                strip.background = element_blank(),
                strip.text.y = element_blank()
            )
    } else {
        stop(cat(sprintf("%s unsupported `strip_position`...\n", strip_position)))
    }

    return(plot)
}


#' Add intron connections to a tibble of exons
#'
#' @param exon_like_tb
#' @param exon_col_start colname of exon/feature starts
#' @param exon_col_end colname of exon/feature ends
#' @param groupby colname of feature group, e.g., "ensembl_transcript_id" (default: NA)
#' This is used to ensure introns only connect between features with same parent
#'
#' @return tibble with jstart, jend added
#' @export
#'
#' @examples
add_introns <- function(
    exon_like_tb,
    exon_col_start="start",
    exon_col_end="end",
    groupby = NA,
    strand_col_name = "strand"
) {
    if (is.na(groupby[1])) {
        groupby <- "rowindex"
        exon_like_tb$rowindex <- 1:nrow(exon_like_tb)
    }
    exon_like_tb_wintrons <-
        exon_like_tb %>%
        dplyr::rename(
            start = exon_col_start,
            end = exon_col_end
        ) %>%
        dplyr::group_by_at(groupby) %>%
        dplyr::arrange(start, end, .by_group=TRUE) %>%
        dplyr::mutate(
            # if positive strand,
            # intron start is exon's end
            # intron end is next exon's start
            # NA if last exon in gene
            # (could be range/feature instead of exon/gene ...)
            jstart = ifelse(
                dplyr::row_number()==dplyr::n(),
                NA,
                end
            ),
            jend = ifelse(
                dplyr::row_number()==dplyr::n(),
                NA,
                dplyr::lead(start, 1)
            )
        ) %>%
        dplyr::ungroup() #%>%

    # tidyr::unite(all_of(groupby), col = "track_group_id")
    return(exon_like_tb_wintrons)
}


#' Add junctions to a plot
#'
#' @param junction_tb <tbl> has jstart, jend, inclorskip, entity
#'  inclorskip: values are 'exclusion' (curve goes down) or 'inclusion' (curves go up)
#' @param gplot genome browser like gplot
#' @param strand <character> "+" or "-" used for defining curvature of a junction
#' (optional)
#' @param junction_colors <character> optional vector of colors where names
#' are correspond to the entity
#' @param exon_height <float> this is how tall the
#' exons are in the plot tracks, so 0.5 takes up half track
#' @param curvature <float> how curved are the junctions (0 is flat)
#'
#' @return
#' @export
#'
#' @examples
plot_junctions <- function(
    junction_tb,
    gplot,
    strand = "+",
    junction_colors = NA,
    exon_height = 0.25,
    curvature = 0.3,
    width = 2
) {
    track_breathing_room_factor <- 2
    gplot <-
        gplot +
        geom_curve(
            data = filter(junction_tb, inclorskip == "exclusion"),
            aes(
                x = jstart,
                y = -exon_height / track_breathing_room_factor,
                xend = jend,
                yend = -exon_height / track_breathing_room_factor,
                # linetype = leafcutter_junction,
                color = entity
            ),
            curvature = curvature * as.integer(paste0(strand, "1")),
            size = width
        ) +
        geom_curve(
            data = filter(junction_tb, inclorskip == "inclusion"),
            aes(
                x = jstart,
                y = exon_height / track_breathing_room_factor,
                xend = jend,
                yend = exon_height / track_breathing_room_factor,
                # linetype = leafcutter_junction,
                color = entity
            ),
            curvature = (-1) * curvature * as.integer(paste0(strand, "1")),
            size = width) +
        theme(legend.title = element_blank()) +
        scale_linetype_discrete(
            guide = F
        )
    if (!is.na(junction_colors[1])) {
        gplot <-
            gplot +
            scale_color_manual(
                guide = F,
                values = junction_colors
            )
    }
    return(gplot)
}

# theme copied from ggpubr source code:
# https://github.com/kassambara/ggpubr/blob/master/R/theme_pubr.R
my_theme_pubclean <- function (base_size = 12, base_family = "", flip = FALSE)
{
  res <- theme_grey(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      legend.position = "top"

    )
  if(flip){
    res <- res + theme(
      panel.grid.major.x = element_line(linetype = "dotted", color = "grey"),
      axis.line.y = element_line(color = "black")
    )
  }
  else{
    res <- res + theme(
      panel.grid.major.y = element_line(linetype = "dotted", color = "grey")
    )
  }
  res
}
