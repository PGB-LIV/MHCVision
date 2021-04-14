# import modules
import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import beta

argv = sys.argv
path = 'support_data/'
# load supported alleles 
supported_allele_netpan = []
supported_allele_flurry = []
with open(path+'supplied_alleles_NetMHCpan.txt', 'rt') as fhla:
    for line in fhla:
        supported_allele_netpan.append(line.strip())
with open(path+'supplied_alleles_MHCflurry.txt', 'rt') as fhla:
    for line in fhla:
        supported_allele_flurry.append(line.strip())
"""
check errors
"""     


# check if the user provided valid arguments
def check_valid_argument(arg):
    invalid_flag = False
    if '-a' not in arg and '--allele' not in arg:
        invalid_flag = True
        print('Error: -a/--allele argument is required')
    if '-t' not in arg and '--tool' not in arg:
        invalid_flag = True
        print('Error: -t/--tool argument is required')
    if '-i' not in arg and '--input' not in arg:
        invalid_flag = True
        print('Error: -i/--input argument is required')
    # print detected error
    if invalid_flag == True:
        print('\nPlease see help information below:')
    return invalid_flag


# print help statement
def print_help_():
    print('usage: mhcvision.py [options] input_file.csv -o/--output output_file.csv\n'
          '-a, --allele   REQUIRED: type the allele name i.e. HLA-A0101, which are supported in the "supplied_alleles.txt"\n'
          '-t, --tool     REQUIRED: Specify the MHC-peptide prediction tool you used, type NetMHCpan or MHCflurry\n'
          '-i, --input    REQUIRED: specify the input filename,\n'
          '               the input file must be in ".CSV" format (comma-separated values)'
          ', the column headers must contain "Peptide", "IC50","%Rank"\n'
          '-o, --output   Optional: specify the output filename\n'
          '-h, --help     Print the usage information')


# extract argument values
def extract_required_arg(arg):
    if '-a' in arg:
        allele_loc = arg.index('-a')
    else:
        allele_loc = arg.index('--allele')
    if '-t' in arg:
        tool_loc = arg.index('-t')
    else:
        tool_loc = arg.index('--tool')
    if '-i' in arg:
        input_loc = arg.index('-i')
    else:
        input_loc = arg.index('--input')
    if '-o' in arg or '--output' in arg:
        out_file = arg[-1]
    else:
        out_file = 'output_' + arg[input_loc+1]
    return arg[allele_loc+1], arg[tool_loc+1], arg[input_loc+1], out_file


# check the input table file format and allele name
def check_input_arg(hla, prediction, file):
    invalid_flag = False
    # 1.check input prediction tools, they must be NetMHCpan or MHCflurry
    if prediction not in ['NetMHCpan', 'MHCflurry']:
        invalid_flag = True
        print('Error: '+ prediction +' has not been supported in this version, it must be NetMHCpan or MHCflurry')
    # 2.check input allele
    supported_allele = []
    if prediction == 'NetMHCpan':
        supported_allele = supported_allele_netpan
    if prediction == 'MHCflurry':
        supported_allele = supported_allele_flurry
    if hla not in supported_allele:
        invalid_flag = True
        print('Error: '+ hla + ' is not in "supplied_alleles_'+prediction+'.txt"')
    # 3.check input file format
    df = pd.read_csv(file, sep=',')
    row1 = list(df.iloc[0,:])
    if len(row1) <= 1:
        invalid_flag = True
        print('Error: The input file must be .CSV format, the columns are separated by ","')
    header = df.columns
    # 4.check the present of IC50 column
    if 'IC50' not in header:
        invalid_flag = True
        print('Error: "IC50" column can not be found')
    # 5.check values in IC50 column
    ic50_col = list(df.loc[:, 'IC50'])
    res = [n for n in ic50_col if (isinstance(n, str))]
    if res:
        invalid_flag = True
        print('Error: "IC50" column must be numbers')
    # print detected error
    if invalid_flag:
        print('\nPlease see help information below:')
    return invalid_flag


"""
The working scripts
"""

# collect values from each iterations
weight = {'k1': [], 'k2': []}
alpha_shape = {'k1': [], 'k2': []}
beta_shape = {'k1': [], 'k2': []}
# initial values
k = 2  # Numbers of components in the mixture
pi = np.zeros(2)
means = np.zeros(2)
variances = np.zeros(2)
a = np.zeros(2)  # alpha
b = np.zeros(2)  # beta


# convert IC50 scores to range of 0-1 for beta distribution
# input file must contain four columns of log IC50 with header;HLA,pep,IC50,%rank
def convert_score_for_beta(file):
    df_template = pd.read_csv(file, index_col=False)
    template = np.log10(df_template.loc[:, 'IC50'])  # take log10
    max_score = np.max(template)
    # convert scores to in range of (0,1)
    converted_template = template / max_score
    for i in range(len(converted_template)):
        if converted_template[i] == 1:
            converted_template[i] = 0.99999  # do not allow score = 1
    return converted_template


class BMM:
    def __init__(self, score, hla, hla_parameter, file):
        # Initialise class method
        self.template = file  # input file containing IC50 scores
        self.hla = hla  # input hla allele
        self.par_range = hla_parameter  # beta parameter ranges
        self.input = score  # data


    """
    1. define the number of cluster and calculate initial parameters
    """
    def initialisation(self):
        input_data = self.input
        size = len(input_data)
        # initial weights
        weights_i = [0.5, 0.5]
        len_r1 = int(round(weights_i[0] * size))
        len_r2 = int(round(weights_i[1] * size))
        d = abs(size - (len_r1 + len_r2))
        if d != 0:
            if len_r1 > len_r2:
                len_r2 += d
            if len_r2 > len_r1:
                len_r1 += d
        # calculate initial parameters
        converted_template = np.sort(input_data)
        # per each component, initial alpha and beta were calculated from mean and variance
        means_i = [np.mean(converted_template[:len_r1]), np.mean(converted_template[len_r1:])]
        variances_i = [np.var(converted_template[:len_r1]), np.var(converted_template[len_r1:])]
        means_i = np.nan_to_num(means_i)
        variances_i = np.nan_to_num(variances_i)
        # update weights, means, variances, a, b
        for z in range(0, k):
            pi[z] = weights_i[z]
            means[z] = means_i[z]
            variances[z] = variances_i[z]
            phi_i = ((means[z] * (1 - means[z])) / variances[z]) - 1
            a_i = phi_i * means[z]
            b_i = phi_i * (1 - means[z])
            a[z] = a_i
            b[z] = b_i
        return converted_template

    """
    2. pdf for mixture beta
    """
    def pdf(self, alp, bet):
        input_data = self.input
        pdf_x = []
        for i in range(0, len(input_data)):
            score = float(input_data[i])
            beta_pdf = beta.pdf(score, alp, bet)
            pdf_x.append(beta_pdf)
        return pdf_x

    """
    3. E steps 
    - do with each component
    - Finds the expected labels of each data point
    """
    def expectation(self):
        likelihood_k = []
        for i in range(0, k):
            mu = means[i]
            va = variances[i]
            phi = ((mu * (1 - mu)) / va) - 1
            al = np.nan_to_num(float(phi * mu))
            be = np.nan_to_num(float(phi * (1 - mu)))
            likelihood = self.pdf(al, be)
            likelihood_k.append(likelihood)
        likelihood_k = np.array(likelihood_k)
        return likelihood_k

    """
    3. M steps 
    - do with each component
    - Finds the maximum likelihood parameters of our model
    - use the current values for the parameters to evaluate the posterior probabilities 
    of the data to have been generated by each distribution
    """
    def maximisation(self, delta):
        hla_parameter_range = self.par_range
        input_data = self.input
        w = []
        likelihood_k = self.expectation()
        mu = means
        va = variances
        wg = pi
        key = ['k1', 'k2']
        # call parameter ranges of an input hla
        a2_b2 = list(hla_parameter_range[self.hla])  # [a2_min,a2_max,b2_min,b2_max]
        min_a2 = a2_b2[0]
        max_a2 = a2_b2[1]
        min_b2 = a2_b2[2]
        max_b2 = a2_b2[3]
        for i in range(0, k):
            wi = np.nan_to_num(likelihood_k[i] * wg[i])  # expected responsibility of each component
            wj = []
            for j in range(0, k):
                wij = likelihood_k[j] * wg[j]  # expected responsibility of two components
                wj.append(wij)
            wj = np.nan_to_num(wj)
            nw = np.nan_to_num(wi / sum(wj))  # expected responsibility of each component i for each data point x, w(ij)
            w.append(nw)
            # update mixture proportions
            pi[i] = np.mean(w[i])  # new mixture proportion --> pi(j)
            # update means and variances
            means[i] = np.sum(w[i] * input_data) / (np.sum(w[i]))
            variances[i] = np.sum(w[i] * np.square(input_data - mu[i])) / (np.sum(w[i]))
            # update alpha and beta from means and variances
            phi = ((mu[i] * (1 - mu[i])) / va[i]) - 1
            a[i] = float(phi * mu[i])
            b[i] = float(phi * (1 - mu[i]))
            # constrain only the second component, a2 and b2
            if i == 1:
                ## when delta < 1e-3, constraining will be start ##
                if delta <= 1e-3:
                    # constrain a and b to be not out of range
                    if a[i] > max_a2:
                        a[i] = max_a2
                    if a[i] == np.nan or a[i] < min_a2:
                        a[i] = min_a2
                    if b[i] > max_b2:
                        b[i] = max_b2
                    if b[i] == np.nan or b[i] < min_b2:
                        b[i] = min_b2
                # update phi, mean, and variance from  a and b
                phi = a[i] + b[i]
                means[i] = a[i] / phi
                variances[i] = (means[i] * (1 - means[i])) / (1 + phi)
            # record data
            weight[key[i]].append(pi[i])
            alpha_shape[key[i]].append(a[i])
            beta_shape[key[i]].append(b[i])
        return

    """
    4. Termination
    - do iterate the EM step
    - each step return the Kq, max relative change of estimate parameters from the previous round
    - stop and return the best estimated parameters when Kq < 1e-5
    """
    def termination(self):
        w1,w2,a1,a2,b1,b2 = 0,0,0,0,0,0
        tq = 1e-5
        kq = 1000
        fix_kq = 10000
        # the fist two iteration
        for n in range(2):
            self.expectation()
            self.maximisation(kq)
        while kq > tq:
            key = ['k1', 'k2']
            c_w1 = list(weight[key[0]])
            c_w2 = list(weight[key[1]])
            c_a1 = list(alpha_shape[key[0]])
            c_a2 = list(alpha_shape[key[1]])
            c_b1 = list(beta_shape[key[0]])
            c_b2 = list(beta_shape[key[1]])
            # find kq of each iteration
            delta_a1 = abs(c_a1[-1] - c_a1[-2]) / max(c_a1[-1], c_a1[-2])
            delta_a2 = abs(c_a2[-1] - c_a2[-2]) / max(c_a2[-1], c_a2[-2])
            delta_b1 = abs(c_b1[-1] - c_b1[-2]) / max(c_b1[-1], c_b1[-2])
            delta_b2 = abs(c_b2[-1] - c_b2[-2]) / max(c_b2[-1], c_b2[-2])
            delta_w1 = abs(c_w1[-1] - c_w1[-2]) / max(c_w1[-1], c_w1[-2])
            delta_w2 = abs(c_w2[-1] - c_w2[-2]) / max(c_w2[-1], c_w2[-2])
            kq = max(delta_a1, delta_b1, delta_w1, delta_a2, delta_b2, delta_w2)
            if kq <= 1e-3:
                fix_kq = 1e-3
            # once got delta less than 1e-3, force the following to go to constraining even their delta > 1e-3
            if kq <= fix_kq:
                self.expectation()
                self.maximisation(kq)
            if kq > fix_kq:
                self.expectation()
                self.maximisation(fix_kq)
            # terminate when kq < tq
            if 0 < kq < tq:
                w1,w2,a1,a2,b1,b2 = c_w1[-1],c_w2[-1],c_a1[-1],c_a2[-1],c_b1[-1],c_b2[-1]
                break
        # write beta_parameter.csv
        with open('beta_parameter_'+self.template, 'wt') as fpar:
            fpar.write('parameter' + ',' + 'value' + '\n' +
                       'pi1' + ',' + str(w1) + '\n' +
                       'pi2' + ',' + str(w2) + '\n' +
                       'a1' + ',' + str(a1) + '\n' +
                       'a2' + ',' + str(a2) + '\n' +
                       'b1' + ',' + str(b1) + '\n' +
                       'b2' + ',' + str(b2) + '\n')
        return


class FDR:
    def __init__(self, parameter, score, hla, int_file, out_file):
        # Initialiser class method
        self.par = parameter  # estimated true and false parameters
        self.input = score  # converted IC50
        self.template = int_file  # input file containing IC50 scores
        self.output = out_file
        self.hla = hla  # input hla allele

    def estimate_fdr_pep(self):
        beta_par = pd.read_csv(self.par)
        w1 = beta_par.iloc[0, 1]
        w2 = beta_par.iloc[1, 1]
        a1 = beta_par.iloc[2, 1]
        a2 = beta_par.iloc[3, 1]
        b1 = beta_par.iloc[4, 1]
        b2 = beta_par.iloc[5, 1]
        scores = sorted(list(self.input))
        size = len(scores)
        size1 = int(round((w1 * size), 0))
        size2 = int(round((w2 * size), 0))
        fdr = []
        pep = []
        inversed_pep = []
        for i in range(0,size):
            s = scores[i]
            # global FDR: cumulative density function of beta distribution
            cdf1 = beta.cdf(s, a1, b1)
            scaled_cdf1 = cdf1 * size1
            cdf2 = beta.cdf(s, a2, b2)
            scaled_cdf2 = cdf2 * size2
            fdr_score = scaled_cdf2 / (scaled_cdf1 + scaled_cdf2)
            fdr.append(fdr_score)
            # local FDR: probability density function of beta distribution
            d1 = beta.pdf(s, a1, b1)
            scaled_d1 = d1 * size1
            d2 = beta.pdf(s, a2, b2)
            scaled_d2 = d2 * size2
            pep_score = scaled_d2 / (scaled_d1 + scaled_d2)
            inv_pep = 1-pep_score  # true posterior probability
            pep.append(pep_score)
            inversed_pep.append(inv_pep)
        return fdr, pep, inversed_pep

    def write_output(self):
        # open real data for being template
        int_file = self.template
        df_template = pd.read_csv(int_file, index_col=False).sort_values(by='IC50')
        fdr, pep, inversed_pep = self.estimate_fdr_pep()
        df_template['FDR'] = fdr
        df_template['PEP'] = pep
        df_template['True probability (1-PEP)'] = inversed_pep
        df_template.to_csv(self.output, index=False)
        return


"""
Run the working script
"""
# check for -h/--help argument
if '-h' in argv or '--help' in argv or check_valid_argument(argv) == True:
    print_help_()

else:
    # if everything has been checked out, extract all arguments
    allele, prediction_tool, input_file, output_file = extract_required_arg(argv)
    # check input argument values
    if check_input_arg(allele, prediction_tool, input_file) == True:
        print_help_()
    # if all argument values are correct, make the estimation
    else:
        # loaded hla parameter ranges based on which prediction tool
        df_parameter_range = pd.read_csv(path+'parameter_range_'+prediction_tool+'.csv')
        hla_list = df_parameter_range.iloc[:,0]
        hla_parameter_range = {}
        for y in range(len(hla_list)):
            a2_range = list(df_parameter_range.iloc[y,1:3]) # [min,max]
            b2_range = list(df_parameter_range.iloc[y,3:])
            values = a2_range+b2_range
            hla_parameter_range[hla_list[y]] = values
        data = convert_score_for_beta(input_file)  # converted IC50
        print('The parameter estimation is running...')
        model = BMM(data, allele, hla_parameter_range, input_file,)
        model.initialisation()
        model.termination()
        est_fdr = FDR('beta_parameter_' + input_file, data, allele, input_file, output_file)
        print('The FDR/PEP are calculating...')
        est_fdr.write_output()
        os.remove('beta_parameter_' + input_file)
        print('Done! Wrote output to ' + output_file)
