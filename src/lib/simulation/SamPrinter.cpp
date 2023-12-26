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



#include "simulation/SamPrinter.hpp"


namespace dragenos {
namespace simulation {

void SamPrinter::printQName(const std::string& seqName, const std::uint64_t refPos, const std::string& cigar )
{
   os_ << seqName << ':' << refPos + 1 << ':' << cigar;

}
std::string SamPrinter::generateQual(const std::string& seq)
{
  std::string ret (seq.length(), 'F');

  return ret;
}

void SamPrinter::printSam(const std::string& refName, std::uint64_t refPos, std::string& cigar, std::string& seq)
{
  static constexpr std::uint32_t flag = 2;
  static constexpr std::uint32_t mapQ = 60;
  static constexpr char rNext = '*';
  static constexpr std::uint32_t pNext = 0;
  static constexpr std::uint32_t tLen = 0;

  printQName(refName, refPos, cigar);
  os_ << '\t' << flag << '\t' << refName << '\t' << refPos + 1 << '\t' << mapQ << '\t' << cigar << '\t' << rNext << '\t'
      << pNext << '\t' << tLen << '\t' << seq << '\t' << generateQual(seq) << '\n';

}


}  // namespace simulation
}  // namespace dragenos


