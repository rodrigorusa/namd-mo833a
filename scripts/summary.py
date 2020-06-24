import os
import csv
import argparse
import pandas as pd

test_cases_name = {
    'TC-1': 'ApoA1',
    'TC-2': 'ATPase',
    'TC-3': 'STMV'
}

configs_name = {
    'CFG-1': 'c5.large-2x',
    'CFG-2': 'c5.large-4x',
    'CFG-3': 'c5.large-8x',
    'CFG-4': 'c5.large-16x'
}

def summary(path):

    # Filter
    lines = [ line for line in open(path + '/namd.out') if '[MO833]' in line]
    data_str = ''.join(lines)
    
    df = pd.DataFrame([x.split(',') for x in data_str.split('\n')])

    if df.shape[0] <= 1:
        return None

    row = []
    total_time_time = '-'
    if not df[(df[0] == '[MO833] real')].empty:
        total_time_time = df[(df[0] == '[MO833] real')].iloc[0][1]
    total_time_main = '-'
    if not df[(df[0] == '[MO833] Total time')].empty:
        total_time_main = df[(df[0] == '[MO833] Total time')].iloc[0][1]
    beta = '-'
    if not df[(df[0] == '[MO833] Beta')].empty:
        beta = df[(df[0] == '[MO833] Beta')].iloc[0][2]
    pi_avg = ''
    if not df[(df[0] == '[MO833] PI avg')].empty:
        pi_avg = df[(df[0] == '[MO833] PI avg')].iloc[0][2]
    n_pi = ''
    if not df[(df[0] == '[MO833] PI avg')].empty:
        n_pi = df[(df[0] == '[MO833] PI avg')].iloc[0][3]

    first = '-'
    if not df[(df[2] == '1')].empty:
        first = df[(df[2] == '1')].iloc[0][3]

    second = '-'
    if not df[(df[2] == '2')].empty:
        second = df[(df[2] == '2')].iloc[0][3]

    # Avg PIs (2,6)
    avg2_6 = 0
    avg2_6_df = df[df[2].isin(['2','3','4','5','6'])]
    if avg2_6_df.shape[0] == 5:
        acc = 0.0
        for i in range(0, 5):
            acc += float(avg2_6_df.iloc[i][3])
        avg2_6 = acc/avg2_6_df.shape[0]

    # Avg PIs (2,11)
    avg2_11 = 0
    avg2_11_df = df[df[2].isin(['2','3','4','5','6','7','8','9','10','11'])]
    if avg2_11_df.shape[0] == 10:
        acc = 0.0
        for i in range(0, 10):
            acc += float(avg2_11_df.iloc[i][3])
        avg2_11 = acc/avg2_11_df.shape[0]

    row.append(total_time_time)
    row.append(total_time_main)
    row.append(beta)
    row.append(pi_avg)
    row.append(n_pi)
    row.append(first)
    row.append(second)
    row.append(round(avg2_6, 6))
    row.append(round(avg2_11, 6))
    
    return row


def main(results_dir):

    summary_data = []
    header = ['test-case', 'cfg', 'data', 'total_time-time', 'total_time-main', 'beta', 'avg-PIs', 'n-PIs', '1st-PI', '2nd-PI', 'avg(2-6)', 'avg(2-11)']
    summary_data.append(header)

    for test_case in os.listdir(results_dir):
        if os.path.isfile(results_dir + '/' + test_case):
            continue
        for cfg in os.listdir(results_dir + '/' + test_case):
            if os.path.isfile(results_dir + '/' + test_case + '/' + cfg):
                continue
            for date in os.listdir(results_dir + '/' + test_case + '/' + cfg):
                row = [test_cases_name[test_case], configs_name[cfg], date]
                path = '/'.join([results_dir, test_case, cfg, date])
                data = summary(path)
                if data is not None:
                    row.extend(data)
                    summary_data.append(row)

    with open(results_dir + '/experimental_results.summary.csv', 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerows(summary_data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarization')
    parser.add_argument('-f', dest='results_dir', type=str, required=True)
    args = parser.parse_args()

    main(args.results_dir)
