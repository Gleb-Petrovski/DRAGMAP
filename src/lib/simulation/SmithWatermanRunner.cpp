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

#include "../../include/simulation/SmithWatermanRunner.hpp"

#include <assert.h>

#include "align/Database.hpp"
#include "align/SmithWaterman.hpp"

namespace dragenos {
namespace simulation {

std::string SmithWatermanRunner::runSW(
    const std::uint8_t*                         queryStart,
    const std::uint8_t*                         queryEnd,
    const reference::HashtableConfig::Sequence& s,
    const std::uint32_t                         start,
    const std::uint32_t                         end,
    std::uint64_t                               refPos,
    std::uint32_t                               tLen) const
{
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  align::SmithWaterman sw(similarityScores, gapInit, gapExtend);
  std::string          result;

  dragenos::align::Database seq;
  referenceDir_.getReferenceSequence().getBases(s.seqStart + start, s.seqStart + end, seq);
  sw.align(queryStart, queryEnd, &seq.front(), &seq.front() + seq.size(), 10, 2, false, result);

  //  referenceDir.getReferenceSequence().getBases(s.seqStart + refPos - startBuffer ,s.seqStart + refPos + readLength + tLen + endBuffer, seq);
  //

  //    std::cerr << convertToString(seq, referenceDir) << std::endl;
  //        std::cerr << convertToString(read, referenceDir) << std::endl;
  //        std::cerr << result << std::endl;
  //        std::cerr << std::endl;

  return result;
}

}  // namespace simulation
}  // namespace dragenos
