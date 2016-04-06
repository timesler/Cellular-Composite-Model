#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _capump_reg();
extern void _h_reg();
extern void _lva_reg();
extern void _nap_reg();
extern void _spike_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," capump.mod");
fprintf(stderr," h.mod");
fprintf(stderr," lva.mod");
fprintf(stderr," nap.mod");
fprintf(stderr," spike.mod");
fprintf(stderr, "\n");
    }
_capump_reg();
_h_reg();
_lva_reg();
_nap_reg();
_spike_reg();
}
