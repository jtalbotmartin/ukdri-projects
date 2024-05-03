file_path <- "analysis/Micro_cbV1_full_cohort_CV_Abeta_cngeneson_pc_mito_Brain_region_Sex_APOE_CD33_group.tsv"
dt_old <- read.delim(file_path)

x <- dt_old %>%
  mutate(de = case_when(padj <= 0.05 & logFC >= 0.25 ~ "Up",
                        padj <= 0.05 & logFC <= -0.25 ~ "Down",
                        TRUE ~ "NS"
  ))
table(x$de)


de_old <- dt_old %>%
  select(c(gene, logFC, pval, padj)) %>%
  mutate(term = "old",
         de = case_when(padj <= 0.05 & abs(logFC) >= 0.25 ~ "de",
                        TRUE ~ "NS"))


file_path <- "analysis/Micro_cbV1_full_cohort_CV_Abeta_cngeneson_pc_mito_Brain_region_Sex_Age_PostMortemInterval_APOE_CD33_group.tsv"
dt_new <- read.delim(file_path)

y <- dt_new %>%
  mutate(de = case_when(padj <= 0.05 & logFC >= 0.25 ~ "Up",
                        padj <= 0.05 & logFC <= -0.25 ~ "Down",
                        TRUE ~ "NS"
  ))
table(y$de)

de_new <- dt_new %>%
  select(c(gene, logFC, pval, padj)) %>%
  mutate(term = "new",
         de = case_when(padj <= 0.05 & abs(logFC) >= 0.25 ~ "de",
                        TRUE ~ "NS")) 

dt <- rbind(de_old, de_new)


label1 <- "old"
label2 <- "new"

dt <- pivot_wider(dt, names_from = term, values_from = c(2:4,6))
dt <- dt %>%
  mutate(de = case_when(de_old == "de" & de_new == "de" ~ "Both",
                        de_old == "de" & de_new != "de" ~ label1,
                        de_old != "de" & de_new == "de" ~ label2,
                        de_old != "de" & de_new != "de" ~ "NS",
  )
  ) %>%
  mutate(de = factor(de, levels = c(label1, label2,
                                    "Both", "NS"))) %>%
  arrange(desc(de))


palette_choice <- paletteer::paletteer_d("ggsci::default_nejm")


ggplot(dt) +
  geom_point(aes(x = logFC_old, y = logFC_new, color = de, alpha = de)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab(label1)+
  ylab(label2)+
  scale_colour_manual(name = NULL,
                      aesthetics = c("colour", "fill"),
                      values = c(palette_choice[1], palette_choice[2], palette_choice[3], "grey"),
                      breaks = c(label1, label2, "Both", "NS")
  ) +
  scale_alpha_manual(values = c(1, 1, 1, 0.01))