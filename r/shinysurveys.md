# Implementing Shinysurveys

## Learning Objectives

In this KB entry, you will:
  - Develop a survey using Shinysurveys
  - Create additional input modes
  - Connect your Shiny Survey to Google Sheets via shinyapps.io

## Setting up Shiny Surveys

This KB entry assumes you have gone through the previous KB entry on R Shiny and will build upon some of these concepts. In order to use `Shinysurveys`, you will need to install and load `Shinysurveys`:

```
install.packages("shinysurveys")
library(shinysurveys)
```

## Developing a survey using Shinysurveys

When using Shinysurveys, the questions are formatted in a structured dataframe. The dataframe has seven slots:
  - `question` - The set of questions you'd like to ask
  - `option` - This can vary with the different input types so we will discuss them as we discuss the various input types
  - `input_type` - This is the type of input (Selection, text response, multiple choice, etc.)
  - `input_id` - A reference ID for this input
  - `dependence` - Which input ID is this question dependent on. If none, then NA
  - `dependence_value` - Which `option` will trigger the dependency
  - `required` - A boolean if this is required to complete the survey

### First example

Let's go ahead and start with a basic example of text response input:

```
questions_df <- data.frame(
  question = "What is your first name?",
  option = "Enter your first name here",
  input_type = "text",
  input_id = "name",
  dependence = NA,
  dependence_value = NA,
  required = TRUE
)
```

Here we have created the dataframe that `shinysurveys` will employ in order to survey a text response question asking for the respondent's name. Now, we need to set-up the rest of the R Shiny script. The `ui` is pretty simple:

```
ui <- fluidPage(
  surveyOutput(questions_df)
)
```

Lots of the formatting is handled by `shinysurveys`. We just need to specify the dataframe that we want to use within the `surveyOutput()` function. The `server` side is a bit more complicated:

```
server <- function(input, output, session) {
  renderSurvey()
  
  observeEvent(input$submit, {
    response_data <- getSurveyData()
    print(response_data)
  })
}
```

Essentailly, we need to render the survey with `renderSurvey()`, then when we click a "Submit" button on our survey, we would like to print that survey response to our R console. The `getSurveyData()` function grabs the survey data and puts it into a data frame.

Lastly, we need to call our Shiny app:

```
shinyApp(ui, server)
```

The app should look like:

<p align="center">
<img src="img/Name_survey.png" width="600">
</p>






 
