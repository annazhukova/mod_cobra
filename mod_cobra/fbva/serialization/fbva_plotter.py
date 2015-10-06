from mod_cobra.fbva.fbva_manager import constraint_reaction_of_interest, get_biomass_dependency
from mod_cobra.fbva.serialization import format_values, format_r_name
from mod_sbml.serialization.plot_manager import initialise_fig, create_subplot, save_fig

__author__ = 'anna'


def get_legend_part(model, r_id, ys):
    r_name = format_r_name(model.reactions.get_by_id(r_id))
    m, M = min(ys), max(ys)
    value = format_values(m, M)
    return "%s %s" % (r_name, value)


def show_flux_variations_same_plot(model, figure_path, r_ids, r_id_var, bm_id, max_bound, log_x_scale=False,
                                   log_y_scale=False, remove_zeros=False, constrained=True, minimized=True):
    xs, r_id2ys = get_biomass_dependency(model, r_id_var, r_ids, bm_id, max_bound, remove_zeros, constrained, minimized)
    legend, data = [], []
    for r_id, ys in r_id2ys.iteritems():
        legend.append(format_r_name(model.reactions.get_by_id(r_id)))
        data.append(ys)
    r_var_name = format_r_name(model.reactions.get_by_id(r_id_var))
    create_plot(xs, data, legend, u'%s flux max bound (micromol/min/gDW)' % r_var_name,
                u'Flux values (micromol/min/gDW)', "", log_x_scale,
                log_y_scale)
    save_fig(figure_path)
    return xs, data, legend


def show_flux_variations_same_plot_separate_bm(model, figure_path, r_id_sets, r_id_var, bm_id, max_bound,
                                               log_x_scale=False,
                                               log_y_scale=False, remove_zeros=True,
                                               constrained=True, minimized=True):
    all_r_ids = {bm_id, r_id_var}
    for r_ids in r_id_sets:
        all_r_ids |= r_ids
    xs, r_id2ys = get_biomass_dependency(model, r_id_var, all_r_ids, bm_id, max_bound, remove_zeros, constrained,
                                         minimized)

    for r_id, ys in sorted(r_id2ys.iteritems(), key=lambda (r_id, ys): [abs(y) for y in ys]):
        logging.info("%s %s" % (ys, format_r_name(model.reactions.get_by_id(r_id))))

    r_var_name = format_r_name(model.reactions.get_by_id(r_id_var))
    num = len(r_id_sets)
    initialise_fig((16, 8 * num), 10)
    i = 1
    for r_ids in r_id_sets:
        legend, data = [], []
        for r_id in r_ids:
            ys = r_id2ys[r_id]
            legend.append(get_legend_part(model, r_id, ys))
            data.append(ys)
        create_subplot(xs, data, legend, "%s flux max bound" % r_var_name, "Fluxes", "Varying %s" % r_var_name,
                       log_x_scale,
                       log_y_scale, plot_num=200 + num * 10 + i)
        i += 2
    save_fig(figure_path)


def show_reaction_variations_same_plot(model, inhibited_rn_id, r_ids, bm_id, reaction_activity_rate, file_path, x_label,
                                       y_label, title, max_bound, log_x_scale=False, log_y_scale=False,
                                       remove_zeros=False,
                                       constrained=False, minimized=False, num_colors=0):
    if inhibited_rn_id:
        model.change_objective([model.reactions.get_by_id(bm_id)])
        model.optimize()
        bound = abs(model.solution.x_dict[inhibited_rn_id] * reaction_activity_rate)
        constraint_reaction_of_interest(model, inhibited_rn_id, bound)
    legend, data, xs = [], [], []
    for r_id in r_ids:
        xs, r_id2ys = get_biomass_dependency(model, r_id, {bm_id}, bm_id, max_bound, remove_zeros, constrained,
                                             minimized)
        ys = r_id2ys[bm_id]
        data.append(ys)
        r = model.reactions.get_by_id(r_id)
        legend.append(r.name.replace(r.id, '').strip())#format_r_name(r))
    initialise_fig((6, 6), 14)
    create_subplot(xs, data, legend, x_label, y_label, title, log_x_scale, log_y_scale, 111, legend_loc='upper left',
                   bb_to_anchor=(0, 1), num_colors=num_colors)
    save_fig(file_path)
    return xs, data, legend


def show_reaction_values(model, var_r_id, r_ids, bm_id, file_path, x_label, y_label, title, max_bound,
                         log_x_scale=False, log_y_scale=False, remove_zeros=False, constrained=False, minimized=False):
    xs, r_id2ys = get_biomass_dependency(model, var_r_id, r_ids, bm_id, max_bound, remove_zeros, constrained,
                                         minimized)
    legend, data = [], []
    for r_id, ys in r_id2ys.iteritems():
        legend.append(format_r_name(model.reactions.get_by_id(r_id)))
        data.append(ys)
    initialise_fig((6, 6), 14)
    create_subplot(xs, data, legend, x_label, y_label, title, log_x_scale, log_y_scale, 111, legend_loc='upper left',
                   bb_to_anchor=(0, 1))
    save_fig(file_path)
    return xs, data, legend
