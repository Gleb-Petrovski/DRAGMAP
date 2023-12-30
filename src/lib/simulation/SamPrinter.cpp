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
#include <iostream>
#include <vector>

#include "reference/ReferenceSequence.hpp"
#include "simulation/SamPrinter.hpp"

namespace dragenos {
namespace simulation {

void SamPrinter::printQName(const std::string& seqName, const std::uint64_t refPos, const std::string& cigar)
{
  os_ << seqName << ':' << refPos + 1 << ':' << cigar;
}
std::string SamPrinter::generateQual(const std::vector<unsigned char>& seq)
{
  std::string ret(seq.size(), 'F');

  return ret;
}
std::string convertToString(const std::vector<unsigned char>& seq)
{
  std::string ret;
  for (const auto& b : seq) {
    ret += reference::ReferenceSequence::decodeBase(b);
  }
  return ret;
}

std::string packCigar(const std::string& cigar)
{
  std::string res;
  int         num = 0;
  for (std::uint32_t i = 1; i < cigar.size(); i++) {
    if (cigar.at(i) == cigar.at(i - 1)) {
      num++;
    } else {
      res += std::to_string(num + 1) + cigar.at(i - 1);
      num = 0;
    }
  }
  res += std::to_string(num + 1) + cigar.at(cigar.size() - 1);
  return res;
}

void SamPrinter::printSam(
    const std::string&                refName,
    const std::uint64_t               refPos,
    const std::string&                flatCigar,
    const std::uint32_t               tLen,
    const std::vector<unsigned char>& seq)
{
  static constexpr std::uint32_t flag  = 2;
  static constexpr std::uint32_t mapQ  = 60;
  static constexpr char          rNext = '*';
  static constexpr std::uint32_t pNext = 0;
  const std::string              cigar = packCigar(flatCigar);

  printQName(refName, refPos, cigar);
  os_ << '\t' << flag << '\t' << refName << '\t' << refPos + 1 << '\t' << mapQ << '\t' << cigar << '\t'
      << rNext << '\t' << pNext << '\t' << tLen << '\t' << convertToString(seq) << '\t' << generateQual(seq)
      << '\n';
}

}  // namespace simulation
}  // namespace dragenos
