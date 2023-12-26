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

#include <fstream>
#include <limits>

#include "include/workflow/SimulateReadsWorkflow.hpp"
#include "workflow/GenHashTableWorkflow.hpp"
#include "workflow/Input2SamWorkflow.hpp"

int main(int argc, char* argv[])
{
  dragenos::options::DragenOsOptions opts;
  dragenos::common::parse_options(argc, argv, opts);
  if (opts.buildHashTable_ || opts.htUncompress_)
  {
    dragenos::common::run(dragenos::workflow::buildHashTable, opts);
  }
  else if(opts.simulateReads_)
  {
    dragenos::common::run(dragenos::workflow::simulateReads, opts);
  }
  else
  {
    dragenos::common::run(dragenos::workflow::input2Sam, opts);
  }
}
