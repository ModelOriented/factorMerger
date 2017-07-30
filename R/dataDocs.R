#' Breast cancer dataset
#'
#' A clinical dataset from The Cancer Genome Atlas Project from
#' the \code{RTCGA.clinical} package (\url{https://github.com/RTCGA/RTCGA.clinical}).
#'
#'
#' The variables are as follows (original column names from the \code{RTCGA.clinical} package are given in brackets):
#'
#' \itemize{
#'   \item \code{time} -- event time / censoring time (originally: times)
#'   \item \code{patientCode} -- patient's id (originally: bcr_patient_barcode)
#'   \item \code{vitalStatus} -- patient status (originally: patient.vital_status)
#'   \item \code{drugName} -- drug name given to the patient (originally: patient.drugs.drug.drug_name)
#'   \item \code{histologicalType} -- patient histological type (originally: patient.stage_event.tnm_categories.pathologic_categories.pathologic_t)
#'   \item \code{pathologicCategory} -- patient pathologic type
#' }
#'
#' @name BRCA
#' @docType data
#' @usage data(BRCA)
#' @keywords data
NULL


#' PISA 2012 dataset
#'
#' A dataset from \url{http://www.oecd.org/pisa/} about student performance in 2012.
#'
#' \itemize{
#'   \item PV1MATH student's performance (Mathematics)
#'   \item PV1READ student's performance (Reading)
#'   \item PV1SCIE student's performance (Science)
#'   \item CNT student's country
#'   \item OCOD1 parent's occupation
#' }
#'
#' @name pisa2012
#' @docType data
#' @usage data(pisa2012)
#' @keywords data
NULL


#' European Social Survey - happiness
#'
#' A dataset from \url{http://www.oecd.org/pisa/} about student performance in 2012.
#'
#' \itemize{
#'   \item happy
#'   \item PV1READ student's performance (Reading)
#'   \item PV1SCIE student's performance (Science)
#'   \item CNT student's country
#'   \item OCOD1 parent's occupation
#' }
#'
#' @name ess
#' @docType data
#' @usage data(ess)
#' @keywords data
NULL
