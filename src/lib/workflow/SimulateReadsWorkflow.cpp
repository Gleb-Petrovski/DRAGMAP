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
#include "sam/SamGenerator.hpp"
#include "simulation/ReadGenerator.hpp"
#include "simulation/ReadPackager.hpp"
#include "simulation/SamPrinter.hpp"
#include "simulation/SmithWatermanValidator.hpp"
#include "simulation/ValidatorThread.hpp"
#include "simulation/VariantGenerator.hpp"
#include "simulation/WorkQueue.hpp"

namespace dragenos {

namespace workflow {
void simulateReads(const dragenos::options::DragenOsOptions& options)
{
  if (!options.simulateReads_) {
    return;
  }
  const reference::ReferenceDir7 referenceDir(
      options.refDir_, options.mmapReference_, options.loadReference_);

  simulation::ReadGenerator    rGen(options.readLength_, options.readSpacing_);
  simulation::VariantGenerator vGen(options.maxVarLen_, options.varSpacingTarget_);
  std::ofstream                os;

  namespace bfs = boost::filesystem;
  if (!options.outputDirectory_.empty()) {
    if (!exists(bfs::path(options.outputDirectory_))) {
      BOOST_THROW_EXCEPTION(common::IoException(
          ENOENT, std::string("Output directory does not exist: ") + options.outputDirectory_));
    }
    if (options.generateSam_) {
      const auto filePath =
          bfs::path(options.outputDirectory_) / (options.outputFilePrefix_ + "Generated" + ".sam");
      os.open(filePath.c_str());
      if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, std::string("Failed to create SAM file: ") + filePath.string() + ": " + strerror(errno)));
      }
      if (options.verbose_) {
        std::cerr << "INFO: writing SAM file to " << filePath << std::endl;
      }
    } else if (options.validate_) {
      const auto filePath = bfs::path(options.outputDirectory_) /
                            (options.outputFilePrefix_ + "SmithWatermanFidelity" + ".csv");
      os.open(filePath.c_str());
      if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, std::string("Failed to create CSV file: ") + filePath.string() + ": " + strerror(errno)));
      }
      if (options.verbose_) {
        std::cerr << "INFO: writing CSV file to " << filePath << std::endl;
      }
    }
  }
  if (options.generateSam_) {
    std::ostream& samFile = os.is_open() ? os : std::cout;
    sam::SamGenerator::generateHeader(
        samFile, referenceDir.getHashtableConfig(), options.getCommandLine(), options.rgid_, options.rgsm_);
    simulation::SamPrinter samOutput(referenceDir, samFile);
    const auto&            seqs = referenceDir.getHashtableConfig().getSequences();
    for (const auto& s : seqs) {
      const auto vars = vGen.generateVariants(0, s.seqLen);
      rGen.generateReads(s, referenceDir, vars, samOutput);
    }
  } else if (options.validate_) {
    std::ostream&                      csvFile = os.is_open() ? os : std::cout;
    simulation::SmithWatermanValidator validator(
        referenceDir, options.startFlank_, options.endFlank_, csvFile);
    simulation::WorkQueue       workQueue;
    simulation::ReadPackager    packager(50, workQueue);
    simulation::ValidatorThread worker(validator, options.readLength_, workQueue);
    std::thread                 thread1(worker);

    const auto& seqs = referenceDir.getHashtableConfig().getSequences();
    for (const auto& s : seqs) {
      const auto vars = vGen.generateVariants(0, s.seqLen);
      rGen.generateReads(s, referenceDir, vars, packager);
    }
    packager.flush();
    thread1.join();
    validator.printResults();
  }
}

}  // namespace workflow
}  // namespace dragenos
