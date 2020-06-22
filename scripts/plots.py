import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

results_dir = '../experimental_results'

test_cases_name = {
    'TC-1': 'ApoA1',
    'TC-2': 'F1atPase',
    'TC-3': 'STMV'
}

configs_name = {
    'CFG-1': 'c5.large-2x',
    'CFG-2': 'c5.large-4x',
    'CFG-3': 'c5.large-8x',
    'CFG-4': 'c5.large-16x'
}

def plot_paramount_time(test_case, path):

    fig = plt.figure(figsize=(35,15))

    for cfg in os.listdir(path):
        if os.path.isfile(path + '/' + cfg):
            continue
        data_y = []
        for date in os.listdir(path + '/' + cfg):
            lines = [ line for line in open(path + '/' + cfg + '/' + date + '/namd.out') if '[MO833]' in line]
            data_str = ''.join(lines)
            
            df = pd.DataFrame([x.split(',') for x in data_str.split('\n')])

            if df.shape[0] <= 1:
                continue

            if df[(df[0] == '[MO833] PI avg')].empty:
                continue

            n_pi = int(df[(df[0] == '[MO833] PI avg')].iloc[0][3])

            y = []
            for i in range(0, n_pi):
                y.append(float(df.iloc[i][3]))

            data_y.append(y)

        if len(data_y) == 0:
            continue

        matrix_y = np.array(data_y)
        mean = matrix_y.mean(0)
        std = matrix_y.std(0)

        data_x = np.arange(1,len(mean)+1)
        plt.errorbar(data_x, mean, yerr=std, label=configs_name[cfg])

    plt.title('Execution time: ' + test_cases_name[test_case], fontsize=16)
    plt.xlim(0,len(data_x))
    plt.xlabel('iteration #', fontsize=16)
    plt.xticks(np.arange(0, len(data_x)+1, 50.0))
    plt.ylabel('execution time (s)', fontsize=16)
    plt.grid(True, ls="-", lw=0.5, color='0.9')
    plt.legend(fontsize='large')
    plt.savefig(path + '/paramount_time.pdf', bbox_inches='tight')
    plt.close()


def plot_relative_performance(test_case, path):

    data_y = []
    for cfg in os.listdir(path):
        if os.path.isfile(path + '/' + cfg):
            continue
        for date in os.listdir(path + '/' + cfg):
            lines = [ line for line in open(path + '/' + cfg + '/' + date + '/namd.out') if '[MO833]' in line]
            data_str = ''.join(lines)
            
            df = pd.DataFrame([x.split(',') for x in data_str.split('\n')])

            if df.shape[0] <= 1:
                continue

            if df[(df[0] == '[MO833] PI avg')].empty:
                continue

            n_pi = int(df[(df[0] == '[MO833] PI avg')].iloc[0][3])

            y = []
            for i in range(0, n_pi):
                y.append(float(df.iloc[i][3]))

            data_y.append(y)

    matrix_y = np.array(data_y)
    min_y = matrix_y.min(0)

    fig = plt.figure(figsize=(35,15))

    for cfg in os.listdir(path):
        if os.path.isfile(path + '/' + cfg):
            continue
        data_y = []
        for date in os.listdir(path + '/' + cfg):
            lines = [ line for line in open(path + '/' + cfg + '/' + date + '/namd.out') if '[MO833]' in line]
            data_str = ''.join(lines)
            
            df = pd.DataFrame([x.split(',') for x in data_str.split('\n')])

            if df[(df[0] == '[MO833] PI avg')].empty:
                continue

            n_pi = int(df[(df[0] == '[MO833] PI avg')].iloc[0][3])

            y = []
            for i in range(0, n_pi):
                y.append(float(df.iloc[i][3])/min_y[i])

            data_y.append(y)

        if len(data_y) == 0:
            continue

        matrix_y = np.array(data_y)
        mean = matrix_y.mean(0)
        std = matrix_y.std(0)

        data_x = np.arange(1,len(mean)+1)
        plt.errorbar(data_x, mean, yerr=std, label=configs_name[cfg])

    plt.title('Relative performance: ' + test_cases_name[test_case], fontsize=16)
    plt.xlim(0,len(data_x))
    plt.xlabel('iteration #', fontsize=16)
    plt.xticks(np.arange(0, len(data_x)+1, 50.0))
    plt.ylabel('relative performance (slowdown)', fontsize=16)
    plt.grid(True, ls="-", lw=0.5, color='0.9')
    plt.legend(fontsize='large')
    plt.savefig(path + '/relative_performance.pdf', bbox_inches='tight')
    plt.close()


def main():

    for test_case in os.listdir(results_dir):
        if os.path.isfile(results_dir + '/' + test_case):
            continue
        plot_paramount_time(test_case, results_dir + '/' + test_case)
        plot_relative_performance(test_case, results_dir + '/' + test_case)

if __name__ == '__main__':
    main()
