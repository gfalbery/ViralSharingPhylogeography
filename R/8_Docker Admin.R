
# Docker ####

library(holepunch)

write_compendium_description(package = "ViralSharingPhylogepgraphy", 
                             description = "Predicting the mammalian viral sharing network")

# to write a description, with dependencies. Be sure to fill in placeholder text

write_dockerfile(maintainer = "gfalbery") 

# To write a Dockerfile. It will automatically pick the date of the last 
# modified file, match it to that version of R and add it here. You can 
# override this by passing r_date to some arbitrary date
# (but one for which a R version exists).

generate_badge() # This generates a badge for your readme.

# ----------------------------------------------
# At this time 🙌 push the code to GitHub 🙌
# ----------------------------------------------

# And click on the badge or use the function below to get the build 
# ready ahead of time.
build_binder()
# 🤞🚀