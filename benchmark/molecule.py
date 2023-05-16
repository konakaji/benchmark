import logging
from qwrapper.hamiltonian import Hamiltonian
from qwrapper.obs import PauliObservable
from openfermion import FermionOperator
from openfermion.chem import MolecularData, geometry_from_pubchem
from openfermion.transforms import get_fermion_operator, bravyi_kitaev, jordan_wigner
from openfermionpyscf import run_pyscf


def parse(operator: FermionOperator, nqubit):
    coeffs = []
    paulis = []
    identity_coeff = 0
    for t in operator:
        for term, coefficient in t.terms.items():
            dict = {}
            for index, p_char in term:
                dict[index] = p_char
                if index > nqubit - 1:
                    raise AttributeError("nqubit is not correct.")
            results = []
            is_identity = False
            if len(dict) == 0:
                is_identity = True
            for q_index in range(nqubit):
                if q_index in dict:
                    results.append(dict[q_index])
                    continue
                results.append("I")
            if not is_identity:
                coeffs.append(coefficient.real)
                p_string = "".join(results)
                paulis.append(PauliObservable(p_string))
            else:
                identity_coeff += coefficient.real
    return coeffs, paulis, identity_coeff


class MolecularHamiltonian(Hamiltonian):
    def __init__(self, nqubit, basis, pubchem_name, bravyi_kitaev=True):
        self.identity_coeff = 0
        self.bravyi_kitaev = bravyi_kitaev
        hamiltonian = self._get_molecule_hamiltonian(basis, pubchem_name)
        coeffs, paulis, identity_coeff = parse(hamiltonian, nqubit)
        self.identity_coeff = identity_coeff
        logging.log(logging.INFO, f"# of terms is {len(coeffs)}")
        super().__init__(coeffs, paulis, nqubit)

    def _get_molecule_hamiltonian(self, basis, pubchem_name):
        geometry = geometry_from_pubchem(pubchem_name)
        multiplicity = 1
        charge = 0
        molecule = MolecularData(geometry, basis, multiplicity, charge)
        logging.info("start running pyscf")
        molecule = run_pyscf(molecule, run_scf=0, run_fci=0)
        logging.info("finish running pyscf.")
        molecule_hamiltonian = molecule.get_molecular_hamiltonian()
        logging.info("start qubit mapping.")
        if self.bravyi_kitaev:
            result = bravyi_kitaev(get_fermion_operator(molecule_hamiltonian))
        else:
            result = jordan_wigner(get_fermion_operator(molecule_hamiltonian))
        logging.info("finish qubit mapping.")
        return result
