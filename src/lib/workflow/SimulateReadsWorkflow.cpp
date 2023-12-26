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

#include "options/DragenOsOptions.hpp"
#include "simulation/ReadGenerator.hpp"
#include "simulation/VariantGenerator.hpp"
#include "sam/SamGenerator.hpp"


namespace dragenos {
namespace workflow {



void simulateReads(const dragenos::options::DragenOsOptions& options)
{
  if (!options.simulateReads_) {
    return;
  }
  const reference::ReferenceDir7 referenceDir(
      options.refDir_, options.mmapReference_, options.loadReference_);

  std::ofstream os;
  namespace bfs = boost::filesystem;
  if (!options.outputDirectory_.empty()) {
    if (!exists(bfs::path(options.outputDirectory_))) {
      BOOST_THROW_EXCEPTION(common::IoException(
          ENOENT, std::string("Output directory does not exist: ") + options.outputDirectory_));
    }
    const auto filePath = bfs::path(options.outputDirectory_) / (options.outputFilePrefix_ + ".sam");
    os.open(filePath.c_str());
    if (!os) {
      BOOST_THROW_EXCEPTION(common::IoException(
          errno, std::string("Failed to create SAM file: ") + filePath.string() + ": " + strerror(errno)));
    }
    if (options.verbose_) {
      std::cerr << "INFO: writing SAM file to " << filePath << std::endl;
    }
  }
  std::ostream& samFile = os.is_open() ? os : std::cout;
  sam::SamGenerator::generateHeader(
      samFile, referenceDir.getHashtableConfig(), options.getCommandLine(), options.rgid_, options.rgsm_);


  simulation::SamPrinter output(samFile);
  simulation::ReadGenerator gen(output);

  const auto& seqs = referenceDir.getHashtableConfig().getSequences();
  const auto& seqNames = referenceDir.getHashtableConfig().getSequenceNames();
  for (const auto& s : seqs){
    const auto vars = simulation::generateVariants(0, s.seqLen);
    gen.generateReads(s, referenceDir, vars, seqNames.at(s.id_));

  }

}

}  // namespace workflow
}  // namespace dragenos


