outputDir <- '/home/everett/i/cnobles-CART19_intsite_manuscript_analysis-1fc2e2c/data'
patient_data <- read.csv(
  file.path(outputDir, "cart19_intsite_sample_list.csv")
) %>%
  dplyr::distinct(Patient_ID, Disease, Response_to_Treatment) %>%
  dplyr::rename(
    "patient" = Patient_ID,
    "disease" = Disease,
    "response" = Response_to_Treatment) %>%
  dplyr::mutate(
    clin_trial = str_extract(patient, "[0-9]+"),
    patient = gsub("_", "-", patient),
    disease = gsub("Adult ALL", "aALL", disease),
    disease = gsub("Pediatric ALL", "pALL", disease),
    disease = factor(disease, levels = c("pALL", "aALL", "CLL")),
    response = gsub(" but relapse", "\nw relapse", response),
    response = gsub(" - with transformed dz", "\nw TnDz", response),
    response = factor(
      response,
      levels = c(
        "None", "Partial", "Partial\nw TnDz",
        "Complete", "Complete\nw relapse")),
    simple_response = factor(
      str_extract(as.character(response), "[\\w]+"),
      levels = c("None", "Partial", "Complete")),
    general_response = factor(ifelse(
      response == "None", "Non-responder", "Responder"),
      levels = c("Non-responder", "Responder")),
    determinant_response = factor(ifelse(
      response %in% c("Complete", "Partial\nw TnDz", "Complete\nw relapse"),
      "CR_PRtd", "PR_NR"), levels = c("CR_PRtd", "PR_NR"))) %>%
  dplyr::arrange(disease, response)
