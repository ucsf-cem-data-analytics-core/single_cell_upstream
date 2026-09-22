#!/usr/bin/env Rscript

# =============================================================================
# 1_copy_cellranger_multi_outputs.R
#
# Copy key Cell Ranger multi outputs (counts, h5, BCR, TCR) from a cellranger
# output directory into an organized run-specific results directory, ready for
# downstream loading into Seurat.
#
# Usage:
#   Rscript 1_copy_cellranger_multi_outputs.R [--dry-run] <cellranger_output_dir> <run_results_dir>
#
# Example:
#   Rscript 1_copy_cellranger_multi_outputs.R \
#     /data/processed_data/cellranger_run_output \
#     /data/processed_data/run_results
#
# =============================================================================

# --- Command-line argument parsing -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# Help / usage message
print_help <- function() {
  cat("
================================================================================
  1_copy_cellranger_multi_outputs.R
================================================================================

  DESCRIPTION:
    Copies key Cell Ranger 'multi' pipeline outputs from a cellranger output
    directory into a clean, organized results directory. This makes it easy
    to load the relevant files (counts matrices, h5 files, VDJ contigs) into
    Seurat or other downstream tools without navigating the deep cellranger
    folder structure.

    The script auto-detects all sample folders inside <cellranger_output_dir>
    (any subdirectory containing an 'outs/' folder) and copies the following
    for each sample:
      - Raw feature barcode matrix files       -> <run_results_dir>/multi_counts/<sample>/
      - Raw & filtered h5 matrix files         -> <run_results_dir>/multi_counts_hd5/<sample>/
      - BCR contig annotations & FASTA         -> <run_results_dir>/multi_bcr/<sample>/
      - TCR contig annotations & FASTA         -> <run_results_dir>/multi_tcr/<sample>/

    Files that already exist in the destination are skipped (safe to re-run).

  USAGE:
    Rscript 1_copy_cellranger_multi_outputs.R [--dry-run] <cellranger_output_dir> <run_results_dir>

  ARGUMENTS:
    cellranger_output_dir   Full path to the Cell Ranger output folder.
                            This is the top-level directory that contains one
                            subfolder per sample, each with an 'outs/' directory
                            (e.g., /data/processed_data/cellranger_run_output).

    run_results_dir         Full path to the results directory where the
                            relevant Cell Ranger files will be copied to.
                            This directory will be created if it does not exist.
                            Organize by run so each experiment has its own
                            results folder for easy Seurat loading
                            (e.g., /data/processed_data/run_results).

  OPTIONS:
    --dry-run, -n           Preview mode. Shows exactly which files would be
                            copied and where, without actually copying anything
                            or creating any directories. Use this first to
                            verify the copy plan looks right before committing.

    --help, -h              Show this help message and exit.

  EXAMPLES:
    # Preview what will be copied (nothing is touched on disk)
    Rscript 1_copy_cellranger_multi_outputs.R --dry-run \\
      /data/processed_data/cellranger_run_output \\
      /data/processed_data/run_results

    # Actually copy
    Rscript 1_copy_cellranger_multi_outputs.R \\
      /data/processed_data/cellranger_run_output \\
      /data/processed_data/run_results

  NOTES:
    - Existing files in the destination are NOT overwritten (safe to re-run).
    - Samples are auto-detected: any subfolder of cellranger_output_dir that
      contains an 'outs/' directory is treated as a sample.
    - The destination subdirectories (multi_counts, multi_counts_hd5,
      multi_bcr, multi_tcr) are created automatically.

================================================================================
\n")
}

# Check for --help / -h flag
if (length(args) > 0 && args[1] %in% c("--help", "-h", "-help")) {
  print_help()
  quit(save = "no", status = 0)
}

# Extract --dry-run / -n flag (can appear anywhere before the positional args)
DRY_RUN <- FALSE
flag_idx <- which(args %in% c("--dry-run", "-n"))
if (length(flag_idx) > 0) {
  DRY_RUN <- TRUE
  args <- args[-flag_idx]
}

# Validate argument count (after removing flags)
if (length(args) != 2) {
  cat("\nERROR: Expected 2 positional arguments but received", length(args), "\n")
  cat("\nUsage: Rscript 1_copy_cellranger_multi_outputs.R [--dry-run] <cellranger_output_dir> <run_results_dir>\n")
  cat("Run with --help for detailed usage information.\n\n")
  quit(save = "no", status = 1)
}

# Parse arguments
cellranger_output_dir <- args[1]
run_results_dir       <- args[2]

# Validate that cellranger_output_dir exists
if (!dir.exists(cellranger_output_dir)) {
  cat("\nERROR: cellranger_output_dir does not exist:\n ")
  cat(cellranger_output_dir, "\n")
  cat("Please provide the full path to your Cell Ranger output folder.\n\n")
  quit(save = "no", status = 1)
}

mode_label <- if (DRY_RUN) "DRY RUN (no files will be copied)" else "LIVE RUN"

cat("\n============================================================\n")
cat("  Cell Ranger Multi Output Copier\n")
cat("  Mode:", mode_label, "\n")
cat("============================================================\n")
cat("  Input (cellranger_output_dir):", cellranger_output_dir, "\n")
cat("  Output (run_results_dir):     ", run_results_dir, "\n")
cat("============================================================\n\n")

# --- Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(stringr)
})

# --- Output directory structure ----------------------------------------------

multi_cellranger_out_dirs <- c("multi_counts","multi_counts_filt","multi_counts_hd5","multi_bcr","multi_tcr")

# --- FUNCTIONS ---------------------------------------------------------------

### generateSubDirectories ----
generateSubDirectories <- function(run_dir, dir_list, dry_run = FALSE) {
  for (dir in dir_list) {
    output_dir <- file.path(run_dir, dir)
    if (dry_run) {
      cat("  [DIR]  ", output_dir, "\n")
    } else {
      if (!file.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    }
  }
}

### copySampleCellrangerOutput ----
copySampleCellrangerOutput <- function(sample, input_dir, output_dir, output.type = "counts", sample_out_name = NULL, dry_run = FALSE) {
  # Get sample results directory path
  out_dir_path <- file.path(input_dir, sample, "outs")

  # Only continue if 'outs' directory exists
  if (file.exists(out_dir_path)) {
    # Detect Cell Ranger version by checking for outs/multi/ directory
    # Pre-v10: outs/multi/count/, outs/multi/vdj_b/, outs/multi/vdj_t/
    # v10+:    outs/ (counts directly), outs/vdj_b/, outs/vdj_t/
    legacy_multi <- file.path(input_dir, sample, "outs", "multi")
    is_v10 <- !dir.exists(legacy_multi)

    if (is_v10) {
      sample_dir_count <- out_dir_path                # outs/
      sample_dir_vdj   <- out_dir_path                # outs/vdj_b, outs/vdj_t
    } else {
      sample_dir_count <- file.path(legacy_multi, "count")  # outs/multi/count/
      sample_dir_vdj   <- legacy_multi                      # outs/multi/vdj_b, outs/multi/vdj_t
    }
    sample_dir_filt  <- file.path(input_dir, sample, "outs", "per_sample_outs", sample)

    target_files <- c()
    target_fhs   <- c()

    if (output.type == "counts") {
      if (file.exists(sample_dir_count)) {
        sample_dir   <- file.path(sample_dir_count, "raw_feature_bc_matrix")
        target_files <- list.files(sample_dir)
        target_fhs   <- paste(sample_dir, target_files, sep = "/")
      }
    }

    if (output.type %in% c("counts_hd5")) {
      if (file.exists(sample_dir_count)) {
        target_file_hd5_raw  <- file.path(sample_dir_count, "raw_feature_bc_matrix.h5")
        target_file_hd5_filt <- file.path(sample_dir_filt, "count", "sample_filtered_feature_bc_matrix.h5")

        target_files <- c("raw_feature_bc_matrix.h5", "sample_filtered_feature_bc_matrix.h5")
        target_fhs   <- c(target_file_hd5_raw, target_file_hd5_filt)
      }
    }

    if (output.type %in% c("bcr", "tcr")) {
      vdj_name <- ""
      if (output.type == "bcr") { vdj_name <- "vdj_b" }
      if (output.type == "tcr") { vdj_name <- "vdj_t" }

      # All contigs
      target_files1 <- c()
      target_fhs1   <- c()
      if (file.exists(file.path(sample_dir_vdj, vdj_name))) {
        sample_dir1   <- file.path(sample_dir_vdj, vdj_name)
        target_files1 <- list.files(sample_dir1)
        target_files1 <- target_files1[grep("_contig_annotations.csv|_contig.fasta$", target_files1)]
        target_fhs1   <- paste(sample_dir1, target_files1, sep = "/")
      }

      # Filtered contigs
      target_files2 <- c()
      target_fhs2   <- c()
      if (file.exists(file.path(sample_dir_filt, vdj_name))) {
        sample_dir2   <- file.path(sample_dir_filt, vdj_name)
        target_files2 <- list.files(sample_dir2)
        target_files2 <- target_files2[grep("_contig_annotations.csv|_contig.fasta$", target_files2)]
        target_fhs2   <- paste(sample_dir2, target_files2, sep = "/")
      }

      # All contigs from outs/ (all_contig.fasta, all_contig_annotations.csv)
      target_files3 <- c()
      target_fhs3   <- c()
      if (file.exists(out_dir_path)) {
        all_contig_candidates <- c("all_contig.fasta", "all_contig_annotations.csv")
        for (ac in all_contig_candidates) {
          ac_path <- file.path(out_dir_path, ac)
          if (file.exists(ac_path)) {
            target_files3 <- c(target_files3, ac)
            target_fhs3   <- c(target_fhs3, ac_path)
          }
        }
      }

      target_files <- c(target_files1, target_files2, target_files3)
      target_fhs   <- c(target_fhs1, target_fhs2, target_fhs3)
    }

    # Determine destination folder
    out_sample_dir      <- output_dir
    sample_out_dir_name <- sample
    if (!is.null(sample_out_name)) {
      sample_out_dir_name <- sample_out_name
    }

    if (output.type == "counts")     { out_sample_dir <- file.path(output_dir, "multi_counts", sample_out_dir_name) }
    if (output.type == "counts_hd5") { out_sample_dir <- file.path(output_dir, "multi_counts_hd5", sample_out_dir_name) }
    if (output.type == "bcr")        { out_sample_dir <- file.path(output_dir, "multi_bcr", sample_out_dir_name) }
    if (output.type == "tcr")        { out_sample_dir <- file.path(output_dir, "multi_tcr", sample_out_dir_name) }

    if (dry_run) {
      # --- DRY RUN: print the copy plan ---
      if (length(target_files) > 0) {
        for (i in 1:length(target_files)) {
          src  <- target_fhs[i]
          dest <- file.path(out_sample_dir, target_files[i])
          if (file.exists(src)) {
            fsize <- file.info(src)$size
            fsize_str <- if (fsize < 1024) {
              paste0(fsize, " B")
            } else if (fsize < 1024^2) {
              paste0(round(fsize / 1024, 1), " KB")
            } else if (fsize < 1024^3) {
              paste0(round(fsize / 1024^2, 1), " MB")
            } else {
              paste0(round(fsize / 1024^3, 2), " GB")
            }
            cat("  [COPY] ", src, "\n")
            cat("      -> ", dest, "  (", fsize_str, ")\n")
          } else {
            cat("  [SKIP] ", src, "  (source file not found)\n")
          }
        }
      } else {
        cat("  [NONE] ", sample, ": no", output.type, "files found\n")
      }
    } else {
      # --- LIVE RUN: actually copy ---
      if (!file.exists(out_sample_dir)) { dir.create(out_sample_dir, recursive = TRUE) }

      if (length(target_files) > 0) {
        for (i in 1:length(target_files)) {
          file      <- target_files[i]
          target_fh <- target_fhs[i]

          if (!file.exists(file.path(out_sample_dir, file))) {
            system2("cp", c(target_fh, out_sample_dir))
          }
        }
      }
    }

  } else {
    if (dry_run) {
      cat("  [SKIP] ", sample, ": no outs/ directory found\n")
    }
  }
}

### copyCellrangerOutput ----
copyCellrangerOutput <- function(samples, input_dir, output_dir, output.type = "counts", dry_run = FALSE) {
  lapply(as.list(samples), copySampleCellrangerOutput,
         input_dir = input_dir, output_dir = output_dir,
         output.type = output.type, dry_run = dry_run)
}

# --- MAIN EXECUTION ----------------------------------------------------------

# Auto-detect samples: any subdirectory containing an "outs/" folder
all_subdirs <- list.dirs(cellranger_output_dir, full.names = FALSE, recursive = FALSE)
samples <- all_subdirs[sapply(all_subdirs, function(s) {
  dir.exists(file.path(cellranger_output_dir, s, "outs"))
})]

if (length(samples) == 0) {
  cat("WARNING: No sample directories found with an 'outs/' subfolder in:\n")
  cat(" ", cellranger_output_dir, "\n")
  cat("Nothing to copy. Exiting.\n\n")
  quit(save = "no", status = 0)
}

cat("Found", length(samples), "sample(s):", paste(samples, collapse = ", "), "\n\n")

if (DRY_RUN) {
  cat("--- Directory structure that would be created ---\n")
  cat("  [DIR]  ", run_results_dir, "\n")
  generateSubDirectories(run_results_dir, multi_cellranger_out_dirs, dry_run = TRUE)
  cat("\n--- Files that would be copied ---\n")
} else {
  # Create output directory structure
  if (!dir.exists(run_results_dir)) {
    dir.create(run_results_dir, recursive = TRUE)
    cat("Created results directory:", run_results_dir, "\n")
  }
  generateSubDirectories(run_results_dir, multi_cellranger_out_dirs)
}

# Copy (or preview) each output type for all samples
output_types <- c("counts", "counts_hd5", "bcr", "tcr")
for (otype in output_types) {
  cat("\n >>", toupper(otype), "\n")
  copyCellrangerOutput(samples,
                       input_dir  = cellranger_output_dir,
                       output_dir = run_results_dir,
                       output.type = otype,
                       dry_run = DRY_RUN)
}

if (DRY_RUN) {
  cat("\n============================================================\n")
  cat("  Dry run complete. No files were copied or created.\n")
  cat("  To execute for real, remove the --dry-run flag.\n")
  cat("============================================================\n\n")
} else {
  cat("\n============================================================\n")
  cat("  Done! Results copied to:", run_results_dir, "\n")
  cat("============================================================\n\n")
}
