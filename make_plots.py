# Usage: python make_plots.py file_path plot_path
import sys
sys.path.append('.')
import hydraulic_jump_2D

if __name__ == "__main__":
    file_path = sys.argv[1]
    plot_path = sys.argv[2]
    hydraulic_jump_2D.plot_surface(file_path,make_anim=False,save_plots=True,plotdir=plot_path)
