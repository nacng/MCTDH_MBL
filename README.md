Folder structure:
|- L
   |- gamma

Data structure:
First row: 1st element is number of averages of the observable. All other entries are irrelevant

For the remaining rows:
Columns are in the following order, where Mean_n is the nth observable (for example, the nth site when calculating Sz; or the occupation of the nth level of the qudit)

Time, Mean_1, Mean_2, ..., Sample std dev_1, Sample std dev_2, ...

where the standard error of the mean is given by 
(Sample std dev_n)/sqrt(N)
with N being the total number of runs (first element of each .csv file)

Observables:
MIavg		Mutual information between two contiguous halves of the spin chain. Plotted using log2, not ln.
savg		Entanglement entropy between half of the spins and the rest of the system. Plotted using log2, not ln.
qdsavg		Entanglement entropy between the qudit and the full spin chain. Computed using log2, not ln.
qdVar		Qudit variance
qdWF		Occupations of the different states of the qudit
sgQ		Spin glass order parameter
sz		Average magnetization of each spin
