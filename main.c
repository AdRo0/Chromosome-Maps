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

    int current_chr = 0;
    gboolean *visited = g_new0(gboolean, genes_count);

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

        // Run mapping for this chromosome
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
            // Single gene: place at 0
            for (int a = 0; a < component->len; a++)
            {
                int i = g_array_index(component, int, a);
                genes[i].pos = 0;
                genes[i].chromosome = current_chr;
            }
            g_array_free(component, TRUE);
            current_chr++;
            continue;
        }

        double d = haldane_to_cm(frequencies[anchor_i][anchor_j]);
        genes[anchor_i].pos = -d / 2.0;
        genes[anchor_j].pos = d / 2.0;
        genes[anchor_i].chromosome = current_chr;
        genes[anchor_j].chromosome = current_chr;

        // Place others
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

            if (dist_i != -1 && dist_j != -1)
            {
                double total = dist_i + dist_j;
                double ratio = dist_i / total;
                genes[k].pos = genes[anchor_i].pos + ratio * (genes[anchor_j].pos - genes[anchor_i].pos);
                dist_i = (fabs(genes[k].pos - genes[anchor_i].pos));
                dist_j = (fabs(genes[k].pos - genes[anchor_j].pos));
                if (haldane_to_r(dist_i) * 100 > frequencies[k][anchor_i] * 1.7 || haldane_to_r(dist_j) * 100 > frequencies[k][anchor_j] * 1.7)
                    return -1;
                else
                {
                    g_print("DIFF I: %.2f / %.2f", haldane_to_r(dist_i), frequencies[k][anchor_i] * 2.0);
                }
            }
            else if (dist_i != -1 && dist_j == -1)
            {
                g_print("DIST I: %s\n", genes[k].name);
                genes[k].pos = genes[anchor_i].pos + dist_i;
                genes[k].possible_pos = genes[anchor_i].pos - dist_i;
            }
            else if (dist_j != -1 && dist_i == -1)
            {
                genes[k].pos = genes[anchor_j].pos - dist_j;
                genes[k].possible_pos = genes[anchor_j].pos + dist_j;
            }
            else
            {
                for (int b = 0; b < genes_count; b++)
                {
                    if (frequencies[k][b] == -1)
                        continue;

                    if (genes[b].pos == -1)
                        continue;

                    genes[k].pos = genes[b].pos + haldane_to_cm(frequencies[k][b]);
                    genes[k].possible_pos = genes[b].pos - haldane_to_cm(frequencies[k][b]);
                    g_print("ELSE: %s, %.2f, %.2f\n", genes[k].name, genes[k].pos, genes[k].possible_pos);
                }
            }

            genes[k].chromosome = current_chr;
        }

        // Resolve ambiguities
        int changes = 1;
        while (changes > 0)
        {
            changes = 0;
            for (int a = 0; a < component->len; a++)
            {
                int i = g_array_index(component, int, a);
                if (genes[i].possible_pos == -1)
                    continue;

                for (int b = 0; b < component->len; b++)
                {
                    int j = g_array_index(component, int, b);
                    if (i == j)
                        continue;
                    if (frequencies[i][j] == -1 || frequencies[i][j] == 50)
                        continue;
                    if (genes[j].pos == -1 || genes[j].possible_pos != -1 || genes[i].possible_pos == -1)
                        continue;
                    if (j == anchor_i || j == anchor_j)
                        continue;

                    double dist = haldane_to_cm(frequencies[i][j]);
                    int match = best_match(genes[j].pos, genes[i].pos, genes[i].possible_pos, dist);
                    g_print("MATCHING\n");
                    if (match == -1)
                    {
                        return -1;
                    }
                    else if (match == 0)
                    {
                        genes[i].possible_pos = -1;
                    }
                    else
                    {
                        genes[i].pos = genes[i].possible_pos;
                        genes[i].possible_pos = -1;
                    }
                    changes++;
                }
            }
        }

        g_array_free(component, TRUE);
        current_chr++;
    }

    g_free(visited);
    return 0;
}

static gboolean on_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data)
{
    double offset_x = pan_x;
    int width = gtk_widget_get_allocated_width(widget);
    int height = gtk_widget_get_allocated_height(widget);

    double rect_w = 900;
    double rect_h = 100;
    double radius = 20.0;

    // Find number of chromosomes
    int max_chr = -1;
    for (int i = 0; i < genes_count; i++)
    {
        if (genes[i].chromosome > max_chr)
            max_chr = genes[i].chromosome;
    }
    int chr_count = max_chr + 1;

    double chr_spacing = 200; // vertical spacing between chromosomes

    for (int chr = 0; chr < chr_count; chr++)
    {
        // Y offset for this chromosome
        double y_offset = chr * chr_spacing + 50;

        // Apply offset_x consistently - only once to the base position
        double x0 = (width - rect_w) / 2.0 + offset_x;
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
            if (genes[i].chromosome != chr)
                continue;
            if (genes[i].pos < min_pos)
                min_pos = genes[i].pos;
            if (genes[i].pos > max_pos)
                max_pos = genes[i].pos;
        }
        if (min_pos == DBL_MAX)
            continue; // empty chromosome

        double span = max_pos - min_pos;
        if (span == 0)
            span = 1;

        double margin = rect_w * 0.05;
        double usable_w = rect_w - 2 * margin;

        // Precompute x positions - these are already relative to x0 which includes offset
        double *gene_x = g_new(double, genes_count);
        for (int i = 0; i < genes_count; i++)
        {
            if (genes[i].chromosome != chr)
                continue;
            double norm = (genes[i].pos - min_pos) / span;
            gene_x[i] = x0 + margin + norm * usable_w;
        }

        // Find most negative gene position to shift positions so left-most = 0
        double min_pos_global = DBL_MAX;
        for (int i = 0; i < genes_count; i++)
        {
            if (genes[i].chromosome != chr)
                continue;
            if (genes[i].pos < min_pos_global)
                min_pos_global = genes[i].pos;
        }
        double shift = (min_pos_global < 0) ? -min_pos_global : 0;

        // Draw gene markers + labels + shifted positions
        for (int i = 0; i < genes_count; i++)
        {
            if (genes[i].chromosome != chr)
                continue;

            // Draw vertical line for gene
            cairo_set_source_rgb(cr, 1, 0, 0);
            cairo_set_line_width(cr, 4.0);
            cairo_move_to(cr, gene_x[i], y0);
            cairo_line_to(cr, gene_x[i], y1);
            cairo_stroke(cr);

            // Draw gene name above chromosome
            cairo_set_source_rgb(cr, 0, 0, 0);
            cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
            cairo_set_font_size(cr, 12.0);
            cairo_text_extents_t extents;
            cairo_text_extents(cr, genes[i].name, &extents);
            double label_x = gene_x[i] - extents.width / 2.0;
            double label_y = y0 - 8; // 8 pixels above
            cairo_move_to(cr, label_x, label_y);
            cairo_show_text(cr, genes[i].name);

            // Draw shifted gene position below chromosome
            char pos_buf[32];
            snprintf(pos_buf, sizeof(pos_buf), "%.2f", genes[i].pos + shift);
            cairo_text_extents(cr, pos_buf, &extents);
            double pos_x = gene_x[i] - extents.width / 2.0;
            double pos_y = y1 + 15; // 15 pixels below
            cairo_move_to(cr, pos_x, pos_y);
            cairo_show_text(cr, pos_buf);
        }

        g_free(gene_x);
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