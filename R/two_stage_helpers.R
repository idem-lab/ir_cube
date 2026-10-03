# Helpers shared by the two-stage scripts (#21). Functions only; sources
# R/validation_scoring.R (functions and settings) for the rules the scoring
# uses: thin_draws(), rho_lookup() and string_seed().

source("R/validation_scoring.R")

# timestamped progress line
report <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush(stdout())
}

# wait until enough memory is available, e.g. before loading a saved fold;
# other jobs share the machine
wait_for_memory <- function(needed_gb = 12, poll_seconds = 180) {
  # off where a scheduler outside R budgets the memory (in a container
  # /proc/meminfo is the host's)
  if (Sys.getenv("IR_CUBE_NO_MEMORY_WAIT") == "1") return(invisible(NA))
  repeat {
    meminfo <- readLines("/proc/meminfo")
    available_gb <- as.numeric(gsub("\\D", "", grep("^MemAvailable", meminfo,
                                                    value = TRUE))) / 1024 ^ 2
    if (available_gb >= needed_gb) return(invisible(available_gb))
    report("%.1f GB available, waiting for %.0f GB", available_gb, needed_gb)
    Sys.sleep(poll_seconds)
  }
}

# peak resident memory of this process so far, from the kernel
peak_memory_gb <- function() {
  status <- readLines("/proc/self/status")
  as.numeric(gsub("\\D", "", grep("^VmHWM", status, value = TRUE))) / 1024 ^ 2
}
