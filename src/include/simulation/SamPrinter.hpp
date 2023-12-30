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
#pragma once
#include <iostream>
#include <vector>

namespace dragenos {
namespace simulation {

class SamPrinter {
public:
  std::ostream& os_;
  SamPrinter(std::ostream& os) : os_(os) {}
  void printSam(
      const std::string&                refName,
      const std::uint64_t               refPos,
      const std::string&                cigar,
      const std::vector<unsigned char>& seq);

private:
  void printQName(const std::string& seqName, const std::uint64_t refPos, const std::string& cigar);
  static std::string generateQual(const std::vector<unsigned char>& seq);
};

}  // namespace simulation
}  // namespace dragenos
