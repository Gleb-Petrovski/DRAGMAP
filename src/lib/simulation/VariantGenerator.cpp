/**
 ** DRAGEN Open Source Software
 ** Copyright (c) 2019-2020 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 **/


#include <assert.h>


#include "simulation/VariantGenerator.hpp"

namespace dragenos {
namespace simulation {



Variants generateVariants(std::uint32_t beginPos, std::uint32_t endPos)
{
  std::uint32_t cPos = beginPos;

  const int varSpacingTarget = 50;
  Variants vars;
  while (cPos < endPos){

    const int randomValue = std::rand();
    Variant v;

    v.refPos_ =  cPos;
    if (randomValue %3 == 0){
      v.refLen_ = 0;
      v.seq_ = "IIIII";
    }
    else if(randomValue %3 == 1)
    {
      v.refLen_ = 1;
      v.seq_ = "X";
    }
    else
    {
      v.refLen_ = 5;
      v.seq_ = "";
    }
    vars.push_back(v);
    cPos += v.refLen_ + randomValue % varSpacingTarget + 1;
  }
  return vars;
}


}  // namespace simulation
}  // namespace dragenos


