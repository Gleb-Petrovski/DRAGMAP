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

#include "align/Database.hpp"
#include "simulation/CigarComparer.hpp"

namespace dragenos {
namespace simulation {

struct Cursor {
  std::uint32_t refPos   = 0;
  std::uint32_t cigarPos = 0;
  std::uint32_t readPos  = 0;
};

std::string unitCigar1 = "SSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";

std::string unitCigar2 = "MMMMMMIIIIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
void        moveCounter(Cursor& a, const std::string& cigar)
{
  char c = cigar.at(a.cigarPos);
  if (c != 'N' && c != 'D') {
    ++a.readPos;
  }
  if (c != 'S' && c != 'I') {
    ++a.refPos;
  }

  ++a.cigarPos;
}
std::uint32_t CigarComparer::countMatches(const std::string& cigar)
{
  std::uint32_t ret = 0;
  for (const char b : cigar) {
    if (b == 'M') {
      ++ret;
    }
  }
  return ret;
}

std::uint32_t CigarComparer::compareCigars(const std::string& cigar1, const std::string& cigar2)
{
  Cursor        a;
  Cursor        b;
  std::uint32_t fidelity = 0;

  while (a.cigarPos != cigar1.length() && b.cigarPos != cigar2.length()) {
    //std::cerr << a.readPos << ':' << a.refPos << ':' <<cigar1.at(a.cigarPos) << '\t' << b.readPos << ':' << b.refPos << ':' << cigar2.at(b.cigarPos) <<std::endl;
    if (cigar1.at(a.cigarPos) != 'M') {
      moveCounter(a, cigar1);
      continue;
    }
    if (cigar2.at(b.cigarPos) != 'M') {
      moveCounter(b, cigar2);
      continue;
    }

    if (a.readPos < b.readPos) {
      moveCounter(a, cigar1);
      continue;
    } else if (b.readPos < a.readPos) {
      moveCounter(b, cigar2);
      continue;
    }

    if (cigar1.at(a.cigarPos) == 'M' && cigar2.at(b.cigarPos) == 'M') {
      fidelity += a.refPos == b.refPos;
    }
    moveCounter(a, cigar1);
    moveCounter(b, cigar2);
  }

  return fidelity;
}

void CigarComparer::unitTest()
{
  std::cerr << compareCigars(unitCigar1, unitCigar2) << std::endl;
}

}  // namespace simulation
}  // namespace dragenos
