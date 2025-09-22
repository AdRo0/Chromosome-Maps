#include <gtk/gtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

GtkWidget *n_spin; // spin for n genes_count
GtkWidget *probabilities_grid;
GtkWidget *main_stack;
GtkWidget *main_window;
GtkWidget *drawing_area;

static double pan_x = 0.0;
static double last_drag_x = 0.0;
static gboolean dragging = FALSE;

typedef struct
{
    char *name;
    double pos;
    double possible_pos; // Other possible position
    int chromosome;
} Gene;
Gene genes[5];

int **frequencies;
int **distances; // in cM
int genes_count = 5;

// A complete map = one possible arrangement of all chromosomes
typedef struct
{
    Gene *placements;     // array of size genes_count (copy of global genes)
    int chromosome_count; // number of chromosomes in this map
} ChromosomeMap;

// Holds all possible solutions (each entry is a ChromosomeMap*)
GPtrArray *all_maps = NULL;

// Free function for ChromosomeMap
static void free_chromosome_map(gpointer data)
{
    ChromosomeMap *map = (ChromosomeMap *)data;
    if (!map)
        return;
    if (map->placements)
        g_free(map->placements);
    g_free(map);
}

// Clone a map deeply (copies placements array)
static ChromosomeMap *clone_map(const ChromosomeMap *src)
{
    ChromosomeMap *copy = g_new0(ChromosomeMap, 1);
    copy->chromosome_count = src->chromosome_count;
    copy->placements = g_new0(Gene, genes_count);
    for (int i = 0; i < genes_count; i++)
    {
        copy->placements[i] = src->placements[i]; // struct copy
        // names are shared (no strdup needed, they live in global `genes`)
    }
    return copy;
}

void setup_genes_arr()
{
    for (int j = 1; j <= genes_count; j++)
    {
        GtkWidget *child = gtk_grid_get_child_at(GTK_GRID(probabilities_grid), j, 0);
        const char *text = gtk_entry_get_text(GTK_ENTRY(child));
        Gene gene;
        // strcpy(gene.name, text);
        gene.name = text;
        gene.pos = -1.0;
        gene.possible_pos = -1.0;
        gene.chromosome = -1;
        genes[j - 1] = gene;
    }
}

// Recombination frequency to map distance in centiMorgan
double haldane_to_cm(int r_percent)
{
    double r = r_percent / 100.0;
    if (r < 0.0 || r > 0.5)
    {
        return -1.0;
    }
    double result = -50.0 * log(1.0 - 2.0 * r);
    return result;
}

// map distance in centiMorgan to Recombination frequency
double haldane_to_r(double cm)
{
    if (cm < 0.0)
    {
        fprintf(stderr, "Error: map distance must be non-negative\n");
        return -1.0;
    }
    return 0.5 * (1.0 - exp(-2.0 * cm / 100.0));
}

void setup_frequencies_table()
{
    frequencies = malloc(genes_count * sizeof(int *));
    for (int i = 0; i < genes_count; i++)
    {
        frequencies[i] = malloc(genes_count * sizeof(int));
    }

    for (int i = 1; i <= genes_count; i++)
    {
        for (int j = 1; j <= genes_count; j++)
        {
            if (i > j)
                continue;
            GtkWidget *child = gtk_grid_get_child_at(GTK_GRID(probabilities_grid), j, i);
            const char *text = gtk_entry_get_text(GTK_ENTRY(child));
            int ivalue;

            if (text[0] == '\0')
            {
                ivalue = -1;
            }
            else
            {
                ivalue = atoi(text);
            }

            if (ivalue > 50 || ivalue < -1)
            {
                g_warning("Aborted draw: Inconsistent probabilities!");
                return;
            }

            frequencies[i - 1][j - 1] = ivalue;
            frequencies[j - 1][i - 1] = ivalue;
        }
    }
}

// Returns 0 if pos_a is better. Return 1 if pos_b is better. Return -1 if none
int best_match(double gene_pos, double pos_a, double pos_b, double desired_dist)
{
    double dist_a = fabs(pos_a - gene_pos);
    double dist_b = fabs(pos_b - gene_pos);

    double diff_a = fabs(dist_a - desired_dist);
    double diff_b = fabs(dist_b - desired_dist);

    g_print("difa: %.2f, difb: %.2f, posa: %.2f, posb: %.2f, post: %.2f, dist: %.2f\n", diff_a, diff_b, pos_a, pos_b, gene_pos, desired_dist);

    if (diff_a <= diff_b && pos_a > -97.0)
    {
        // if (diff_a > 10.0)
        //     return -1;
        return 0;
    }
    else if (pos_b < 97.0)
    {
        // if (diff_b > 10.0)
        //     return -1;
        return 1;
    }
    else
        return -1;
}

void create_probs_table()
{
    char buf[16];
    for (int i = 0; i <= genes_count; i++)
    {
        for (int j = 0; j <= genes_count; j++)
        {
            GtkWidget *to_attach;
            if (i == 0 && j == 0)
                continue;
            if (i == 0 && j != 0)
            {
                to_attach = gtk_entry_new();
                gtk_entry_set_width_chars(GTK_ENTRY(to_attach), 6);
                gtk_entry_set_max_length(GTK_ENTRY(to_attach), 5);
                snprintf(buf, sizeof(buf), "G%d", j);
                gtk_entry_set_text(GTK_ENTRY(to_attach), buf);

                g_object_set_data(G_OBJECT(to_attach), "node-index", GINT_TO_POINTER(j));
                gtk_grid_attach(GTK_GRID(probabilities_grid), to_attach, j, i, 1, 1);
                continue;
            }
            if (j == 0 && i != 0)
            {
                to_attach = gtk_entry_new();
                gtk_entry_set_width_chars(GTK_ENTRY(to_attach), 6);
                gtk_entry_set_max_length(GTK_ENTRY(to_attach), 5);
                snprintf(buf, sizeof(buf), "G%d", i);
                gtk_entry_set_text(GTK_ENTRY(to_attach), buf);

                GtkWidget *col_entry = gtk_grid_get_child_at(GTK_GRID(probabilities_grid), i, 0);

                // g_signal_connect(to_attach, "changed", G_CALLBACK(on_name_changed), col_entry);
                // g_signal_connect(col_entry, "changed", G_CALLBACK(on_name_changed), to_attach);

                gtk_grid_attach(GTK_GRID(probabilities_grid), to_attach, j, i, 1, 1);
                continue;
            }
            to_attach = gtk_entry_new();
            gtk_entry_set_width_chars(GTK_ENTRY(to_attach), 6);
            gtk_entry_set_max_length(GTK_ENTRY(to_attach), 5);
            gtk_entry_set_input_purpose(GTK_ENTRY(to_attach), GTK_INPUT_PURPOSE_DIGITS);

            if (i == j)
            {
                gtk_entry_set_text(GTK_ENTRY(to_attach), "0");
                gtk_widget_set_sensitive(to_attach, FALSE);
            }
            else if (i > j)
            {
                gtk_entry_set_text(GTK_ENTRY(to_attach), "-");
                gtk_widget_set_sensitive(to_attach, FALSE);
            }
            else
                gtk_entry_set_text(GTK_ENTRY(to_attach), "");
            gtk_grid_attach(GTK_GRID(probabilities_grid), to_attach, j, i, 1, 1);
        }
    }
    gtk_widget_show_all(probabilities_grid);
}

int calculate_maps(void)
{
    setup_genes_arr();
    setup_frequencies_table();

    if (all_maps)
    {
        g_ptr_array_free(all_maps, TRUE); // free old results if any
    }
    all_maps = g_ptr_array_new_with_free_func(free_chromosome_map);

    // ---- Create the base map (empty placements) ----
    ChromosomeMap *base = g_new0(ChromosomeMap, 1);
    base->placements = g_new0(Gene, genes_count);
    for (int i = 0; i < genes_count; i++)
    {
        base->placements[i].name = genes[i].name; // reuse name pointer
        base->placements[i].pos = -1;
        base->placements[i].possible_pos = -1;
        base->placements[i].chromosome = -1;
    }
    base->chromosome_count = 0;
    g_ptr_array_add(all_maps, base);

    gboolean *visited = g_new0(gboolean, genes_count);
    int current_chr = 0;

    for (int start = 0; start < genes_count; start++)
    {
        if (visited[start])
            continue;

        // Find connected component
        GQueue queue = G_QUEUE_INIT;
        g_queue_push_tail(&queue, GINT_TO_POINTER(start));
        GArray *component = g_array_new(FALSE, FALSE, sizeof(int));

        while (!g_queue_is_empty(&queue))
        {
            int idx = GPOINTER_TO_INT(g_queue_pop_head(&queue));
            if (visited[idx])
                continue;
            visited[idx] = TRUE;
            g_array_append_val(component, idx);

            for (int j = 0; j < genes_count; j++)
            {
                if (j == idx || visited[j])
                    continue;
                if (frequencies[idx][j] != -1 && frequencies[idx][j] < 50)
                {
                    g_queue_push_tail(&queue, GINT_TO_POINTER(j));
                }
            }
        }

        if (component->len == 0)
        {
            g_array_free(component, TRUE);
            continue;
        }

        // For each current map in all_maps, expand it with placements for this chromosome
        int maps_before = all_maps->len;
        for (int m = 0; m < maps_before; m++)
        {
            ChromosomeMap *map = g_ptr_array_index(all_maps, m);

            // Find anchor
            int anchor_i = -1, anchor_j = -1;
            int max_freq = -1;
            for (int a = 0; a < component->len; a++)
            {
                int i = g_array_index(component, int, a);
                for (int b = a + 1; b < component->len; b++)
                {
                    int j = g_array_index(component, int, b);
                    if (frequencies[i][j] > max_freq && frequencies[i][j] < 50)
                    {
                        max_freq = frequencies[i][j];
                        anchor_i = i;
                        anchor_j = j;
                    }
                }
            }

            if (anchor_i == -1 || anchor_j == -1)
            {
                // Single gene: just place at 0
                for (int a = 0; a < component->len; a++)
                {
                    int i = g_array_index(component, int, a);
                    map->placements[i].pos = 0;
                    map->placements[i].chromosome = current_chr;
                }
                map->chromosome_count++;
                continue;
            }

            double d = haldane_to_cm(frequencies[anchor_i][anchor_j]);
            map->placements[anchor_i].pos = -d / 2.0;
            map->placements[anchor_j].pos = d / 2.0;
            map->placements[anchor_i].chromosome = current_chr;
            map->placements[anchor_j].chromosome = current_chr;

            // Place others (branch on ambiguities)
            GPtrArray *new_maps = g_ptr_array_new_with_free_func(free_chromosome_map);
            for (int a = 0; a < component->len; a++)
            {
                int k = g_array_index(component, int, a);
                if (k == anchor_i || k == anchor_j)
                    continue;

                double dist_i = (frequencies[anchor_i][k] != -1 && frequencies[anchor_i][k] < 50)
                                    ? haldane_to_cm(frequencies[anchor_i][k])
                                    : -1;
                double dist_j = (frequencies[anchor_j][k] != -1 && frequencies[anchor_j][k] < 50)
                                    ? haldane_to_cm(frequencies[anchor_j][k])
                                    : -1;

                if (dist_i != -1 && dist_j == -1)
                {
                    // Two possible positions -> branch
                    for (int mi = 0; mi < all_maps->len; mi++)
                    {
                        ChromosomeMap *cur = g_ptr_array_index(all_maps, mi);
                        // copy branch
                        ChromosomeMap *branch = clone_map(cur);
                        branch->placements[k].pos = cur->placements[anchor_i].pos - dist_i;
                        branch->placements[k].chromosome = current_chr;
                        g_ptr_array_add(new_maps, branch);

                        // original map keeps +dist
                        cur->placements[k].pos = cur->placements[anchor_i].pos + dist_i;
                        cur->placements[k].chromosome = current_chr;
                    }
                }
                else if (dist_j != -1 && dist_i == -1)
                {
                    for (int mi = 0; mi < all_maps->len; mi++)
                    {
                        ChromosomeMap *cur = g_ptr_array_index(all_maps, mi);
                        ChromosomeMap *branch = clone_map(cur);
                        branch->placements[k].pos = cur->placements[anchor_j].pos + dist_j;
                        branch->placements[k].chromosome = current_chr;
                        g_ptr_array_add(new_maps, branch);

                        cur->placements[k].pos = cur->placements[anchor_j].pos - dist_j;
                        cur->placements[k].chromosome = current_chr;
                    }
                }
                else
                {
                    // place mid or skip for simplicity
                    for (int mi = 0; mi < all_maps->len; mi++)
                    {
                        ChromosomeMap *cur = g_ptr_array_index(all_maps, mi);
                        cur->placements[k].pos = cur->placements[anchor_i].pos + dist_i;
                        cur->placements[k].chromosome = current_chr;
                    }
                }
            }

            // Merge in new branched maps
            for (int mi = 0; mi < new_maps->len; mi++)
            {
                g_ptr_array_add(all_maps, g_ptr_array_index(new_maps, mi));
            }
            g_ptr_array_free(new_maps, FALSE);

            map->chromosome_count++;
        }

        g_array_free(component, TRUE);
        current_chr++;
    }

    g_free(visited);
    return 0;
}

static gboolean on_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data)
{
    if (!all_maps || all_maps->len == 0)
        return FALSE;

    double offset_x = pan_x;
    int width = gtk_widget_get_allocated_width(widget);
    int height = gtk_widget_get_allocated_height(widget);

    double rect_w = 900;
    double rect_h = 100;
    double radius = 20.0;
    double chr_spacing = 200;          // vertical spacing between chromosomes
    double map_spacing = rect_w + 150; // horizontal spacing between maps

    int map_count = all_maps->len;

    for (int m = 0; m < map_count; m++)
    {
        ChromosomeMap *map = g_ptr_array_index(all_maps, m);

        // Horizontal offset for this map
        double map_x_offset = m * map_spacing;

        // Find number of chromosomes in this map
        int max_chr = -1;
        for (int i = 0; i < genes_count; i++)
        {
            if (map->placements[i].chromosome > max_chr)
                max_chr = map->placements[i].chromosome;
        }
        int chr_count = max_chr + 1;

        int line_index = 0;
        for (int chr = 0; chr < chr_count; chr++)
        {
            // Y offset for this chromosome
            double y_offset = chr * chr_spacing + 50;

            // Base rectangle position
            double x0 = (width - rect_w) / 2.0 + offset_x + map_x_offset;
            double y0 = y_offset;
            double x1 = x0 + rect_w;
            double y1 = y0 + rect_h;

            // Draw chromosome rectangle
            cairo_new_sub_path(cr);
            cairo_arc(cr, x1 - radius, y0 + radius, radius, -90 * (M_PI / 180), 0);
            cairo_arc(cr, x1 - radius, y1 - radius, radius, 0, 90 * (M_PI / 180));
            cairo_arc(cr, x0 + radius, y1 - radius, radius, 90 * (M_PI / 180), 180 * (M_PI / 180));
            cairo_arc(cr, x0 + radius, y0 + radius, radius, 180 * (M_PI / 180), 270 * (M_PI / 180));
            cairo_close_path(cr);

            cairo_set_source_rgb(cr, 0.9, 0.9, 0.9);
            cairo_fill_preserve(cr);
            cairo_set_source_rgb(cr, 0, 0, 0);
            cairo_set_line_width(cr, 2.0);
            cairo_stroke(cr);

            // Collect min/max positions for scaling
            double min_pos = DBL_MAX, max_pos = -DBL_MAX;
            for (int i = 0; i < genes_count; i++)
            {
                if (map->placements[i].chromosome != chr)
                    continue;
                if (map->placements[i].pos < min_pos)
                    min_pos = map->placements[i].pos;
                if (map->placements[i].pos > max_pos)
                    max_pos = map->placements[i].pos;
            }
            if (min_pos == DBL_MAX)
                continue; // empty chromosome

            double span = max_pos - min_pos;
            if (span == 0)
                span = 1;

            double margin = rect_w * 0.05;
            double usable_w = rect_w - 2 * margin;
            double gy = (y0 + y1) / 2.0;

            // Precompute x positions for all genes on this chromosome
            double *gene_x = g_new0(double, genes_count);
            for (int i = 0; i < genes_count; i++)
            {
                if (map->placements[i].chromosome != chr)
                    continue;
                double norm = (map->placements[i].pos - min_pos) / span;
                gene_x[i] = x0 + margin + norm * usable_w;
            }

            // Draw gene markers + labels
            for (int i = 0; i < genes_count; i++)
            {
                if (map->placements[i].chromosome != chr)
                    continue;

                // Red marker line
                cairo_set_source_rgb(cr, 1, 0, 0);
                cairo_set_line_width(cr, 4.0);
                cairo_move_to(cr, gene_x[i], y0);
                cairo_line_to(cr, gene_x[i], y1);
                cairo_stroke(cr);

                // Gene name
                cairo_set_source_rgb(cr, 0, 0, 0);
                cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
                cairo_set_font_size(cr, 12.0);
                cairo_text_extents_t extents;
                cairo_text_extents(cr, map->placements[i].name, &extents);
                double label_x = gene_x[i] - extents.width / 2.0;
                double label_y = y0 - 8;
                cairo_move_to(cr, label_x, label_y);
                cairo_show_text(cr, map->placements[i].name);
            }

            // Draw cM distances for pairs on same chromosome
            for (int i = 0; i < genes_count; i++)
            {
                if (map->placements[i].chromosome != chr)
                    continue;
                for (int j = i + 1; j < genes_count; j++)
                {
                    if (map->placements[j].chromosome != chr)
                        continue;
                    if (frequencies[i][j] == -1 || frequencies[i][j] == 50)
                        continue;

                    double cm = haldane_to_cm(frequencies[i][j]);
                    double offset = (line_index + 1) * 20.0;
                    double yline = y1 + offset;

                    cairo_set_source_rgb(cr, 0, 0, 0);
                    cairo_set_line_width(cr, 1.5);
                    cairo_move_to(cr, gene_x[i], yline);
                    cairo_line_to(cr, gene_x[j], yline);
                    cairo_stroke(cr);

                    char buf[32];
                    snprintf(buf, sizeof(buf), "%.1f cM", cm);
                    double midx = (gene_x[i] + gene_x[j]) / 2.0;
                    cairo_move_to(cr, midx - 10, yline - 4);
                    cairo_show_text(cr, buf);

                    line_index++;
                }
            }

            g_free(gene_x);
        }
    }

    return FALSE;
}

void on_continueBtn_clicked(GtkButton *button, gpointer user_data)
{
    genes_count = (int)gtk_spin_button_get_value(GTK_SPIN_BUTTON(n_spin));
    create_probs_table();
    gtk_stack_set_visible_child_name(GTK_STACK(main_stack), "page1");
}

void on_buildBtn_clicked(GtkButton *button, gpointer user_data)
{
    int maps_result = calculate_maps();
    if (maps_result == -1)
    {
        g_print("FAILED: CHECK ENTRIES!");
        return;
    }
    for (int i = 0; i < genes_count; i++)
    {
        g_print("%s: %2.f / %2.f, ", genes[i].name, genes[i].pos, genes[i].possible_pos);
    }
    gtk_widget_queue_draw(drawing_area);
    gtk_stack_set_visible_child_name(GTK_STACK(main_stack), "page2");
}

static gboolean on_button_press(GtkWidget *widget, GdkEventButton *event, gpointer user_data)
{
    if (event->button == GDK_BUTTON_PRIMARY)
    {
        dragging = TRUE;
        last_drag_x = event->x;
        return TRUE;
    }
    return FALSE;
}

static gboolean on_button_release(GtkWidget *widget, GdkEventButton *event, gpointer user_data)
{
    if (event->button == GDK_BUTTON_PRIMARY)
    {
        dragging = FALSE;
        return TRUE;
    }
    return FALSE;
}

static gboolean on_motion_notify(GtkWidget *widget, GdkEventMotion *event, gpointer user_data)
{
    if (dragging)
    {
        double dx = event->x - last_drag_x;
        pan_x += dx;
        last_drag_x = event->x;
        gtk_widget_queue_draw(widget);
        return TRUE;
    }
    return FALSE;
}

int main(int argc, char *argv[])
{
    gtk_init(&argc, &argv);

    GtkBuilder *builder = gtk_builder_new_from_file("main.glade"); // Si se abre desde el menu

    main_window = GTK_WIDGET(gtk_builder_get_object(builder, "hWindow"));

    gtk_builder_connect_signals(builder, NULL);

    // exit
    g_signal_connect(main_window, "destroy", G_CALLBACK(gtk_main_quit), NULL);

    main_stack = GTK_WIDGET(gtk_builder_get_object(builder, "mainStack"));
    n_spin = GTK_WIDGET(gtk_builder_get_object(builder, "n_spin"));
    probabilities_grid = GTK_WIDGET(gtk_builder_get_object(builder, "probabilities_grid"));
    drawing_area = GTK_WIDGET(gtk_builder_get_object(builder, "drawing_area"));
    g_signal_connect(drawing_area, "draw", G_CALLBACK(on_draw_event), NULL);

    gtk_widget_add_events(drawing_area,
                          GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK | GDK_POINTER_MOTION_MASK);

    g_signal_connect(drawing_area, "draw", G_CALLBACK(on_draw_event), NULL);
    g_signal_connect(drawing_area, "button-press-event", G_CALLBACK(on_button_press), NULL);
    g_signal_connect(drawing_area, "button-release-event", G_CALLBACK(on_button_release), NULL);
    g_signal_connect(drawing_area, "motion-notify-event", G_CALLBACK(on_motion_notify), NULL);

    gtk_widget_show_all(main_window);

    gtk_main();

    g_object_unref(builder);

    return 0;
}