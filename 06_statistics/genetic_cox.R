library("survival")
library("survminer")

data <- read.csv("/home/srd6051/recent_parity_quest.csv")

genes = c(
    "brca1",
    "brca2",
    "palb2",
    "chek2",
    "atm",
    "tp53",
    "pten",
    "stk11",
    "cdh1"
)

data <- data[!is.na(data$years_since_pregnancy),]

# Initialize empty df to record results
results <- data.frame(
    row.names = toupper(genes), 
    HR = rep(0, length(genes)), 
    p.value = rep(0, length(genes))
)

# years_since_pregnancy: Number of years between the last pregnancy and breast cancer diagnosis
# age_at_most_recent_pregnancy: Patient age at the time of the last pregnancy i.e. age at t=0 for Cox regression
# fam_hx: Binary variable as to whether patient has family history or not

for (gene in genes) {
    # Dynamically define formula with gene name
    cox_formula <- as.formula(paste("Surv(years_since_pregnancy) ~", gene, "+ age_at_most_recent_pregnancy + fam_hx"))
    # Fit Cox regression
    res.cox <- coxph(cox_formula, data = data)
    # Obtain and record
    coefs = summary(res.cox)$coefficients
    results[toupper(gene), "HR"] = coefs[gene,"exp(coef)"]
    results[toupper(gene), "p.value"] = coefs[gene,"Pr(>|z|)"]
}

