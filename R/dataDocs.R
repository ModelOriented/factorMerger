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
#' Data from the \code{PISA2012lite} package.
#'
#' The original dataset (Program for International Student Assesment, \url{http://www.oecd.org/pisa/}) contains information on students'
#' performance on various cognitive tests expressed with so-called plausible values measured in three fields: Mathematics (\code{PV1MATH}),
#' Science (\code{PV1SCIE}) and Reading (\code{PV1READ}).
#'
#' \code{pisa2012} is a reweighted version of the \code{PISA2012lite} data. The weighting was performed using the \code{W_FSTUWT} variable with the seed 123.
#'
#' The variables are as follows (original column names from the \code{PISA2012lite} package are given in brackets):
#'
#' \itemize{
#'   \item \code{math} -- student's performance in Mathematics (originally: PV1MATH)
#'   \item \code{reading} -- student's performance in Reading (originally: PV1READ)
#'   \item \code{science} -- student's performance in Science (originally: PV1SCIE)
#'   \item \code{country} -- student's country (originally: CNT)
#' }
#'
#' @name pisa2012
#' @aliases pisaEuro
#' @docType data
#' @usage data(pisa2012)
#' @keywords data
NULL


#' European Social Survey - happiness
#'
#' A dataset from \url{http://www.europeansocialsurvey.org/download.html?file=ESS7e02_1&y=2014} on hapiness
#' of European countries.
#'
#' \code{ess} is a reweighted version of the original data. The weighting was performed using the \code{dweight} variable with the seed 123.
#' The variables are as follows:
#' \itemize{
#'   \item \code{happy} -- binary, denoting happiness of an invidual. Interviewees were asked
#'   to assess their happiness in the scale 0-10 (variable \code{happy}). Respondents who gave value
#'   greater than 5 are considered as "happy",
#'   \item \code{country} -- respondent's country.
#' }
#'
#' @name ess
#' @docType data
#' @usage data(ess)
#' @keywords data
NULL
