import json
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator, FuncFormatter

class PredictData:
    '''Prediction Data for plotting a single individual.

    This class contains the basic information needed to plot the risk of 
    terminal event for a given individual.

    Attributes:
        rec_times: list of floats with recurrent event times.
        pred_point: float with the time where the prediction is made.
        pred_times: list of floats with the end of each prediction interval.
        pred: list of floats with predicted risk of terminal event.
        CI_lower: list of floats with lower end of the Confidence Interval.
        CI_lower: list of floats with upper end of the Confidence Interval.
        display_name: str with legend name for the individual.
    '''

    # Constants
    DEF_DISPLAY_NAMES = ['A',
                         'B',
                         'C',
                         'D',
                         'E']

    def __init__(self, rec_times, pred_point, pred_times, pred, CI_lower,
                 CI_upper, display_name):
        
        # Standarize when there is no recurrent events
        if rec_times:
            self.rec_times = rec_times
        else:
            self.rec_times = []

        self.pred_point = pred_point
        self.pred_times = pred_times
        self.pred = pred
        self.CI_lower = CI_lower
        self.CI_upper = CI_upper
        self.display_name = display_name

    @classmethod
    def from_dict(cls, data_dict):
        '''Turn dictionary with keys into object.
        
        Args:
            data_dict: dictionary with keys as in __init__ arguments.

        Returns:
            A PredictData object.

        '''

        return PredictData(**data_dict)


def shared_frailty_plot(infile, outfile, frmt_options={}):
    '''Plot the risk of terminal event predicted by a shared frailty model.

    Args:
        infile: str with path with JSON file containing plotting information.
        outfile: str with path to export the plot.
        frmt_options: dict with key-value pairs specifying formatting options.
        Valid keys are the ones in DEFAULT_FRMT.

    Returns:
        True if the plot was successfully saved.
    '''

    # Default formatting arguments
    DEFAULT_FRMT = {'time_label': 'time',
                    'rec_label': 'Recurrent event number',
                    'pred_label': 'Terminal event risk between $t$ and $t + w$',
                    'x_minor_scale': None,
                    'x_major_scale': None,
                    'text_plot': True,
                    'pred_ticksize': 'x-small',
                    'lgnd_fontsize': 'small',
                    'max_lgnd_cols': 4,
                    'sc_marker': 'x',
                    'pred_lim': None,
                    'ci_alpha': 0.25,
                    'dpi': 800
    }
    
    # Join the two dictionaries (second overwrites first)
    FRMT = DEFAULT_FRMT | frmt_options 

    # Extra margins
    REC_VER_MAR = 0.05 # For recurrent counter Y-axis
    TEXT_HOR_MAR = -0.015 #-0.5 # For text (horizontal)
    TEXT_VER_MAR = -0.1 # For text (vertical)
                
    # Load data
    with open(infile) as f:
        dat = json.load(f)

    cases = [PredictData.from_dict(case) for case in dat]

    # Info for plotting
    MAX_REC = max(max([len(case.rec_times) for case in cases]), 1) 
    MAX_TIME = max([max(case.pred_times) for case in cases])
    PRED_POINT = cases[0].pred_point

    # Start figure with two Y-axis
    fig, ax_rec = plt.subplots()
    ax_pred = ax_rec.twinx()

    # Plot recurrent event times, predicted risk, and CI for each observation
    rec_sc = []
    pred_c = []
    ci_f = []
    nodname_count = 0
    for case in cases:

        # Asign defaults if there is no display name
        if not case.display_name:
            case.display_name = PredictData.DEF_DISPLAY_NAMES[nodname_count]
            nodname_count += 1

        # Recurrent event
        rec_sc.append(ax_rec.scatter(case.rec_times, 
                                    list(range(1, len(case.rec_times) + 1)),
                                    marker = FRMT['sc_marker'])
        )

        # Prediction (plot returns a list)
        pred_c = pred_c + ax_pred.plot(case.pred_times, case.pred, 
                                       label = case.display_name)
        

        # Confidence interval (same color as line)
        if case.CI_lower and case.CI_upper:
            ci_f.append(ax_pred.fill_between(case.pred_times, case.CI_lower,
                                             case.CI_upper,
                                             color=pred_c[-1].get_color())
                        )

    # Plot prediction point line
    pred_l = ax_rec.vlines(PRED_POINT, ymin=0, ymax= MAX_REC + REC_VER_MAR)

    # Descriptive text
    if FRMT['text_plot']:
        txt = ax_rec.text(PRED_POINT + TEXT_HOR_MAR * MAX_TIME,
                    TEXT_VER_MAR * MAX_REC,
                    r'$t | \longleftrightarrow | t + w$',
                    usetex = True
                    )

    # Formatting

    # X-axes
    ax_rec.set_xlabel(FRMT['time_label'], loc = 'right', usetex = True)
    ax_rec.set_xlim((0, MAX_TIME))

    # Major ticks and labels
    if FRMT['x_major_scale']:
        xtick_frmt = FuncFormatter(lambda x, pos: 
                                   '{:.0f}'.format(x/FRMT['x_major_scale']))

        ax_rec.xaxis.set_major_locator(MultipleLocator(FRMT['x_major_scale']))
        ax_rec.xaxis.set_major_formatter(xtick_frmt)

    # Minor ticks
    if FRMT['x_minor_scale']:    
        ax_rec.xaxis.set_minor_locator(MultipleLocator(FRMT['x_minor_scale']))
    

    # Recurrent axes
    ax_rec.set_ylabel(FRMT['rec_label'], usetex = True)
    ax_rec.set_ylim(0, MAX_REC + REC_VER_MAR)
    ax_rec.set_yticks(list(range(0, MAX_REC + 1)))

    # Prediction axes
    ax_pred.set_ylabel(FRMT['pred_label'], usetex = True)
    ax_pred.set_ylim(bottom=0)
    if FRMT['pred_lim']: ax_pred.set_ylim(top=FRMT['pred_lim'])
    ax_pred.tick_params(axis = 'y', labelsize = FRMT['pred_ticksize'])

    # CI region
    plt.setp(ci_f, 'alpha', DEFAULT_FRMT['ci_alpha'])

    # Prediction line
    pred_l.set_color('black')
    pred_l.set_linestyle(':')

    # Legend
    ax_pred.legend(bbox_to_anchor = (0.5, 1.02),  # Above the plot box
                   loc = 'lower center',
                   ncol = min(len(cases), FRMT['max_lgnd_cols']),
                   fontsize = FRMT['lgnd_fontsize'])
    
    #Save and close
    fig.savefig(outfile, dpi=FRMT['dpi'])
    plt.close(fig)

    return True