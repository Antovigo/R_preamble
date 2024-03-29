line_colour <- 'grey40'
axis_line_weight <- 0.3
grid_line_weight <- 0.1
gridline <- 'grey60'

text_colour <- 'black'
text_size <- 12

half_line <- text_size / 2
spacing <- 1

my_theme <- theme(
	line = element_line(
		colour   = line_colour,
		linewidth = grid_line_weight,
		linetype = 1,
		lineend  = "butt"
	),

	rect = element_rect(
		fill     = "white",
		colour   = text_colour,
		linewidth = 0.5,
		linetype = 1
	),

	text = element_text(
		family     = 'sans',
		face       = "plain",
		colour     = text_colour,
		size       = text_size,
		lineheight = 0.9,
		hjust      = 0.5,
		vjust      = 0.5,
		angle      = 0,
		margin     = margin(),
		debug      = FALSE
	),

	axis.line         = element_line(colour = line_colour, linewidth = axis_line_weight),
	axis.line.x       = NULL,
	axis.line.y       = NULL,
	axis.text         = element_text(size   = rel(0.8), colour = text_colour),
	axis.text.x       = element_text(margin = margin(t = 0.8 * spacing * half_line / 2), vjust = 1),
	axis.text.x.top   = element_text(margin = margin(b = 0.8 * spacing * half_line / 2), vjust = 0),
	axis.text.y       = element_text(margin = margin(r = 0.8 * spacing * half_line / 2), hjust = 1),
	axis.text.y.right = element_text(margin = margin(l = 0.8 * spacing * half_line / 2), hjust = 0),
	axis.ticks        = element_line(colour = line_colour),
	axis.ticks.length = unit(half_line * spacing / 2, "pt"),

	axis.title.x = element_text(
		margin = margin(t = half_line*spacing),
		vjust = 1
	),
	axis.title.x.top = element_text(
		margin = margin(b = half_line*spacing),
		vjust = 0
	),
	axis.title.y = element_text(
		angle = 90,
		margin = margin(r = half_line*spacing),
		vjust = 1
	),
	axis.title.y.right = element_text(
		angle = -90,
		margin = margin(l = half_line*spacing),
		vjust = 0
	),

	legend.background     = element_rect(colour = 'white', fill = 'white'),
	legend.spacing        = unit(0.4*spacing, "cm"),
	legend.spacing.x      = NULL,
	legend.spacing.y      = NULL,
	legend.margin         = margin(0, 0, 0, 0, "cm"),
	legend.key            = element_rect(colour = 'white', fill = 'white'),
	legend.key.size       = grid::unit(1.2, "lines"),
	legend.key.height     = NULL,
	legend.key.width      = NULL,
	legend.text           = element_text(size = rel(0.8)),
	legend.text.align     = NULL, 
	legend.title          = element_text(hjust = 0),
	legend.title.align    = NULL,
	legend.position       = "right",
	legend.direction      = NULL,
	legend.justification  = "center",
	legend.box            = NULL,
	legend.box.margin     = margin(0, 0, 0, 0, "cm"),
	legend.box.background = element_blank(),
	legend.box.spacing    = unit(0.4*spacing, "cm"),

	panel.background      = element_rect(fill = 'white', colour = NA),
	panel.border          = element_blank(),
	panel.grid            = element_line(colour = "white"),
	panel.grid.major      = element_line(linetype = 'dashed', colour = gridline),
	panel.grid.minor      = element_blank(), 
	panel.spacing         = unit(half_line*spacing, "pt"),
	panel.spacing.x       = NULL,
	panel.spacing.y       = NULL,
	panel.ontop           = FALSE,

	strip.background = element_rect(fill = 'white', colour = NA),

	strip.text = element_text(
		colour = text_colour,
		size   = rel(0.8),
		margin = margin(half_line*spacing, half_line*spacing, half_line*spacing, half_line*spacing)
	),

	strip.text.x = NULL,
	strip.text.y = element_text(angle = -90),

	strip.placement   =  "inside",
	strip.placement.x =  NULL,
	strip.placement.y =  NULL,

	strip.switch.pad.grid = unit(0.1, "cm"),
	strip.switch.pad.wrap = unit(0.1, "cm"),

	plot.background = element_rect(
		colour = 'white',
		fill   = 'white'
	),

	plot.title = element_text(
		size   = rel(1.2),
		margin = margin(b = half_line * 1.2 * spacing),
		hjust  = 0,
		vjust  = spacing,
		face   = 'bold'
	),

	plot.subtitle = element_text(
		size   = rel(0.9),
		margin = margin(b = half_line * 0.9 * spacing),
		hjust  = 0,
		vjust  = spacing
	),

	plot.caption = element_text(
		size   = rel(0.9),
		margin = margin(t = half_line * 0.9 * spacing),
		hjust  = spacing,
		vjust  = spacing
	),

	plot.tag = element_text(
		size = rel(1.2),
		hjust = 0.5,
		vjust = 0.5
	),

	plot.tag.position =  'topleft',

	plot.margin   = margin(half_line*spacing, half_line*spacing, half_line*spacing, half_line*spacing),

	complete      = TRUE
)

theme_set(my_theme)
