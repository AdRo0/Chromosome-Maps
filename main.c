#include <gtk/gtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

GtkWidget *n_spin; // spin for n genes
GtkWidget *probabilities_grid;
GtkWidget *main_stack;
GtkWidget *main_window;

int genes = 5;

void create_probs_table()
{
    char buf[16];
    for (int i = 0; i <= genes; i++)
    {
        for (int j = 0; j <= genes; j++)
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
                gtk_entry_set_text(GTK_ENTRY(to_attach), "0");
            gtk_grid_attach(GTK_GRID(probabilities_grid), to_attach, j, i, 1, 1);
        }
    }
    gtk_widget_show_all(probabilities_grid);
}

void on_continueBtn_clicked(GtkButton *button, gpointer user_data)
{
    create_probs_table();
    gtk_stack_set_visible_child_name(GTK_STACK(main_stack), "page1");
}

void on_buildBtn_clicked(GtkButton *button, gpointer user_data)
{
    gtk_stack_set_visible_child_name(GTK_STACK(main_stack), "page1");
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

    gtk_widget_show_all(main_window);

    gtk_main();

    g_object_unref(builder);

    return 0;
}