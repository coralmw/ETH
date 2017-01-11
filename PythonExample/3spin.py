import qutip
import qutip.operators as op
import numpy as np
import matplotlib.pyplot as plt

sigmas = [op.sigmax(), op.sigmay(), op.sigmaz()]

# create the fully connected hamiltonion
Hconnected = op.qzero([2, 2, 2])
for sigma in sigmas:
    for connection in [(0, 1), (1, 2), (2, 0)]:
        Hconnected += qutip.tensor([0.5*sigma if i in connection else op.identity(2) for i in range(3)])
    # for site in range(3):
    #     Hconnected += qutip.tensor([0.5*sigma if i == site else op.identity(2) for i in range(3)])

#Hconnected /= 4.

# we also need the bath hamiltionion - for the 1+2, this is to get singlet/triplet states
Hbath = op.qzero([2, 2])
for sigma in sigmas:
    Hbath += qutip.tensor(0.5*sigma, 0.5*sigma)
    # for site in range(2):
    #     Hbath += qutip.tensor([0.5*sigma if i == site else op.identity(2) for i in range(2)])

bathstates = Hbath.eigenstates()

# system up, bath |ud + du>
psi = qutip.tensor(qutip.basis(2,0), bathstates[1][2])

# calculate the density matrix at each timepoint and then take partial trace
# do for lots of timepoints so we can compare to symbolic results
t = np.linspace(0,2*np.pi, 100)
redens = np.zeros((len(t), 2, 2))
for i, ti in enumerate(t):
    psi_t = sum([ psi.overlap(state) * np.exp( 1j * energy * ti ) * state for energy, state in zip(*Hconnected.eigenstates()) ])
    density_t = psi_t * psi_t.dag()
    redens[i] = density_t.ptrace(0).data.toarray()

# plot it all
redens_overt = redens.transpose((1, 2, 0)) # 2x2's over time
for spin1 in range(2):
    for spin2 in range(2):
        s1name = 'up' if spin1 == 0 else 'down'
        s2name = 'up' if spin2 == 0 else 'down'
        plt.plot(redens_overt[spin1][spin2], label="spin 1 state: {} spin 2 state: {}".format(s1name, s2name))

def aupup(t):
    return 1./9. * (5. + 4.*np.cos(3.*t/2.))

def adowndown(t):
    return 8./9. * (np.sin(3.*t/4.) ** 2)

plt.plot(aupup(t), linestyle='dashed', label="correct upup state")
plt.plot(adowndown(t), linestyle='dashed', label="correct downdown state") # identical
# but the off-diagonal terms are zero

plt.legend()
plt.show()
