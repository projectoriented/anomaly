#!/usr/bin/env bats

# executed before each test
setup() {
    DIR="$( cd "$( dirname "$BATS_TEST_FILENAME" )" >/dev/null && pwd )"
    PATH="$DIR/../workflow/scripts:$PATH"
}

@test "test slurm immediate submit" {
  slurm_immediate_submit.sh 1234567 1234568
}