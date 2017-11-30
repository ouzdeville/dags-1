/*
Implementation by the Keccak, Keyak and Ketje Teams, namely, Guido Bertoni,
Joan Daemen, Michaël Peeters, Gilles Van Assche and Ronny Van Keer, hereby
denoted as "the implementer".

For more information, feedback or questions, please refer to our websites:
http://keccak.noekeon.org/
http://keyak.noekeon.org/
http://ketje.noekeon.org/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef _KetjeMn_h_
#define _KetjeMn_h_

#include "../kecaak/Ketje-common.h"

#ifndef KeccakP800_excluded
    #include "../keccak/KeccakP-800-SnP.h"
    KCP_DeclareKetjeStructure(KetjeMn, KeccakP800_stateSizeInBytes, KeccakP800_stateAlignment)
    KCP_DeclareKetjeFunctions(KetjeMn)
#endif

#endif
