import numpy as np


def burn_detector(multichain, chain_i, chain_j, traces=["posterior", "likelihood", "prior"]):
    # Calculates the burnin that is common between the given chains i and j
    # This is the maximum for i and j (could be indiviual but for convenience we take a common burnin)
    is_log = [tr in multichain[chain_i].log_data for tr in traces]
    if not all(is_log):
        raise ValueError(f"Given trace indicator -{traces[is_log.index(False)]}- is not present in the logfile!")

    burn_ins = []  # storing the burnin for each parameter
    index75 = int(len(
        multichain[chain_i].log_data[traces[0]].values) * 0.75)  # calculating the 75% index, assuming both chains, and all log-parameters have same length
    M = 10  # smoothing parameter

    for tr in traces:
        cur_burn_i, cur_burn_j = index75, index75
        mean_i = np.mean(multichain[chain_i].log_data[tr].values[index75:])
        std_i = np.std(multichain[chain_i].log_data[tr].values[index75:])
        mean_j = np.mean(multichain[chain_j].log_data[tr].values[index75:])
        std_j = np.std(multichain[chain_j].log_data[tr].values[index75:])
        i_set, j_set = False, False
        for i in range(M, index75):
            if not i_set:
                if mean_i - std_i <= np.mean(multichain[chain_i].log_data[tr].values[i-M:i]) <= mean_i + std_i:
                    cur_burn_i = i
                    i_set = True
            if not j_set:
                if mean_j - std_j <= np.mean(multichain[chain_j].log_data[tr].values[i-M:i]) <= mean_j + std_j:
                    cur_burn_j = i
                    j_set = True
            if i_set and j_set:
                break
        burn_ins.append(np.max((cur_burn_j, cur_burn_i)))
    return np.max(burn_ins)
