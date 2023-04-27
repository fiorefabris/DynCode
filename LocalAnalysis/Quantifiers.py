import numpy as np
from itertools import permutations


def Permutation_Entropy(Time_Series, n, Tau, **kwargs):
    """Calculate the Permutation Entropy,  y otras cantidades tambi�n!

    Args:
        Time_Series: Time series for analysis
        n: Order of permutation entropy
        Tau: Time delay between samples, 1 is the classic one. Check [1] eq (4) for reference
    Returns:
        Vector containing Permutation Entropy"""

    N = len(Time_Series)
    Permutations = np.array(list(permutations(range(n))))
    BP_Freq = [0] * len(Permutations)

    for i in range(N - Tau * (n - 1)):
        Ranked_Window_Index = np.array(np.argsort(Time_Series[i:i + Tau * n:Tau], kind='quicksort'))

        Done = False
        j = 0
        while j < len(Permutations) and Done == False:
            if abs(Permutations[j] - Ranked_Window_Index).any() == False:
                BP_Freq[j] = BP_Freq[j] + 1
                Done = True
            j = j + 1

    BP_P = np.divide(np.array(BP_Freq), float(sum(BP_Freq)))

    Aux_BP = [element for element in BP_P if element != 0]
    PE = -sum(Aux_BP * np.log2(Aux_BP)) / np.log2(np.math.factorial(n))

    F = []
    for prev, elem in zip(BP_P[:-1], BP_P[1:]):
        F.append(np.square(np.sqrt(elem) - np.sqrt(prev)))
    F = sum(F)

    P_e = [1 / len(Permutations)] * len(Permutations)
    Aux = np.divide(BP_P + P_e, 2)
    J = -sum(Aux * np.log2(Aux)) + -sum(Aux_BP * np.log2(Aux_BP)) / 2 + -sum(P_e * np.log2(P_e)) / 2
    Q0 = (-2) * 1 / ((N + 1) / N * np.log2(N + 1) - np.log2(2 * N) + np.log2(N))
    C = PE * Q0 * J

    return (BP_P, PE, F, C,Permutations)


def Data_Permutation_Entropy(Time_Series, N, Step):
    """Calculate the Permutation Entropy of a Data Time Serie
Args:
    Time_Series    : Time series for analysis
    N              : Lenght of permutation entropy (list)
    Step           : Time delay between samples, 1 is the classic one.
Returns:
    List containing Permutation Entropy for each value of N
    On the borders, PE is zero"""

    List_PE,List_F,List_C,List_BPP = [],[],[],[]


    for k in N:
        PE,F,C,BPP = [],[],[],[0]*6


        for i in range(len(Time_Series) - Step * (k - 1)):
            Time_Series_W = Time_Series[i:i + k:Step]
            Aux0, Aux1, Aux2, Aux3,Permutations = Permutation_Entropy(Time_Series_W, 3, 1, norm=True)

            BPP = np.add(np.asarray(Aux0), BPP)

            PE.append(Aux1);F.append(Aux2);C.append(Aux3)

        PE = [None] * int(k / 2) + PE + [None] * int((k + 1) / 2 - 1)
        F = [None] * int(k / 2) + F + [None] * int((k + 1) / 2 - 1)
        C = [None] * int(k / 2) + C + [None] * int((k + 1) / 2 - 1)

        List_BPP.extend(BPP); List_PE.extend(PE);List_F.extend(F);List_C.extend(C);

    return (List_BPP, List_PE, List_F, List_C,Permutations)


#%%

def Calculate_FanoF(Function, N):
    n = len(Function)

    Var, Mean = [0] * n,[0] * n

    # PARA EL MEDIO
    for i in range(N // 2, n - N // 2):
        Time_Series = Function[i - N // 2:i + N // 2]
        Var[i] = np.var(Time_Series)
        Mean[i] = np.mean(Time_Series)

    # BORDE IZQ1UIERDO
    for i in range(0, N // 2):
        Time_Series = Function[0:i + N // 2]
        Var[i] = np.var(Time_Series)
        Mean[i] = np.mean(Time_Series)

    # BORDE DERECHO
    for i in range(n - N // 2, n):
        Time_Series = Function[i - N // 2:n]
        Var[i] = np.var(Time_Series)
        Mean[i] = np.mean(Time_Series)

    return (Var, Mean)

#%%
def Calculate_Quantifiers(conc,traze_value,local_time_window_len,delay):
    name = traze_value +'_N'+str(local_time_window_len)
    for C in conc:
        df = conc[C]
        df[name + '_LOCAL_VAR'], df[name + '_LOCAL_MEAN'] = Calculate_FanoF(df[traze_value], local_time_window_len)
        df[name + '_CV2'] = [x**2/y**2 for (x,y) in zip(df[name + '_LOCAL_VAR'],df[name + '_LOCAL_MEAN'])]

        auxPE,auxF,auxC = [],[],[]


        for cells, data in df.groupby(level='cell'):
            Time_Series = data[traze_value]
            List_BPP, List_PE, List_F, List_C,Permutations  = Data_Permutation_Entropy(Time_Series, [local_time_window_len], delay)

            auxPE.extend(List_PE);auxF.extend(List_F);auxC.extend(List_C)

        df[name + '_PE'] = np.array(auxPE); df[name + '_PE'] = 1-df[name + '_PE'].fillna(1)
        df[name + '_F']  = np.array(auxF)
        df[name + '_C']  = np.array(auxC)

    return(Permutations,conc)