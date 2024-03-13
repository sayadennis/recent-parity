library("survival")
library("survminer")

data <- read.csv("/home/srd6051/recent_parity_quest.csv")

genes = c(
    "any_patho_mutation",
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

for (gene in genes) {
    cox_formula <- as.formula(paste("Surv(years_since_pregnancy) ~", gene, "+ age_at_most_recent_pregnancy + fam_hx"))
    res.cox <- coxph(cox_formula, data = data)
    print("")
    print(toupper(gene))
    print(res.cox)
    print("")
}

