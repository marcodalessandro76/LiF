"""
Figures of the article (PumpProbe-NL-Chi/main.tex), made from the runs in NL_Chi, which are on the cluster only.

The data and the calculations are NOT duplicated here: they come from the figure script of the technical notes
(Notes/make_figures.py, loaded as the module `notes`), so that the article and the notes show the same numbers. This
script only sets the style of the article (column width, fonts) and the layout of the panels. A new calculation goes in
Notes/make_figures.py (as a data function, e.g. `hsex_ac_data`, `ip_maps_data`) and is plotted here.

Usage, on ismhpc from ~/work/LiF/PumpProbe-NL-Chi:

    LIF_NL_CHI=$HOME/work/LiF/NL_Chi ~/miniconda3/bin/python make_figures.py [figure ...] [--out DIR]

Without figure names all the figures are made. As for the notes, write the jpeg files outside the cluster copy of the
repository (e.g. --out ~/tmp_claude/paper_figures) and copy them to PumpProbe-NL-Chi/Figures on the laptop, where they
are committed.
"""
import os
import argparse
import importlib.util
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def load_notes():
    # the figure script of the notes as the module `notes` (not with `import make_figures`: this file has the same name)
    spec = importlib.util.spec_from_file_location('notes', os.path.join(HERE, '..', 'Notes', 'make_figures.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


notes = load_notes()

# Style of the article (Phys. Rev. B, revtex4-2 reprint: column width 3.375 in, page width 7.0 in), set after loading
# the notes, which set their own rcParams
COLUMN, PAGE = 3.375, 7.0
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7,
                     'ytick.labelsize': 7, 'savefig.dpi': 300, 'savefig.bbox': 'tight'})

# name: function(out), filled as the figures of the article are defined (Outline.md, "Figures")
FIGURES = {}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Figures of the article on LiF')
    parser.add_argument('figures', nargs='*', help='names of the figures (default: all): %s' % ', '.join(FIGURES))
    parser.add_argument('--out', default=os.path.join(HERE, 'Figures'), help='output folder of the jpeg files')
    args = parser.parse_args()
    out = os.path.expanduser(args.out)
    os.makedirs(out, exist_ok=True)
    for name in args.figures or FIGURES:
        FIGURES[name](out)
