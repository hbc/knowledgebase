# GLOBUS

HMS and FASRC both provide the use of Globus to share/receive data to/from HMS o2 and FASRC. (and even between the two!)  (in globus lingo, HMS_O2/FASRC are considered server-class machines such as campus computing clusters and lab servers)

Depending on what you’re doing, there may be a few steps for setting up a globus share for a client and then giving them instructions. 

## Copying data with a link

However, if you’re just copying data using globus, (for instance from Harvard BPF core), just click on the provided link and follow the Globus instructions (login, specify where you’re copying the data – on HMS RC or FAS RC (ie. /n/data1/cores/bcbio/PIs/PIname/project_folder/data)), click Start> under the data being shared and the transfer should start. 

## More detailed instructions

To have someone send you data via Globus, first create a globus ID here  [https://www.globusid.org/](https://www.globusid.org/). The sender will need this to estalish a transfer.  It will look something like this. e1187b4cb7ce18b2@harvard.edu

When the transfer is ready, you will receive an email message from Globus with a url link to your data.

Create a download destination directory, for example, on o2 in /n/data1/cores/bcbio/PIs/

Click on the globus url link provided in the email message.

Highlight the files/directories in the left hand pane you want to transfer from the client

Click on "Transfer or Sync to"

This will open two panes

In the right hand destination panel, search for "HMS-RC" in the top Collection slot.

When it finds it, paste the path to your dir you've created for the destination directory. Click on the window to enter the path name.

(You may instead have to click on the square to the left of the base directory, like "/n" and click the arrow to progressively open the available subdirectories to your destination directory.)

You may have to select the provided data folders again in the left pane.

You may get a message asking you to confirm your globus id, and enter your Harvard passkey information.

Hit Start, you'll see a window drop down,  and in a moment, a green flag on the left Activity button, as well as a Green "Transfer request submitted successfully" banner, indicating things are working. If you get a red flag something is not right. Stop the transfer and try again.

You can click on the green banner to monitor the transfer.


## Sharing data with a client
Create a Globus guest collection in order to share/receive data from a client

To create a Globus guest collection, navigate to the Globus File Manager, select the collection containing the data you want to share, and click "Share" then "Add Guest Collection". You'll need to give the guest collection a display name and optionally a description. Finally, set permissions for users or groups to access the guest collection. Detailed instruction:

1. [Log into Globus](https://app.globus.org/) and navigate to the [File Manager](https://app.globus.org/file-manager) (top tab on menu to the left).
2. Select the collection that has the files/folders you wish to share (/n/data1/cores/bcbio/PIs/PIname/project_folder/data).
3. Choose the folder that you would like to share and Click Share in the right command pane.
3a.If Share is not available, - you can only share a folder not a file - sometimes you’ll need to create a specific “share” folder and put the file/files there
4. Click the Add Guest Collection button.
5. Provide a name for the guest collection (20250606-PIname project hbc number - anything particular to data) and click Create Collection. (You may need to authenticate again to share/add permissions - use your HMS O2 login)
6. When your collection is created, you’ll be taken to the Permissions tab, where you can set permissions. 
6a.Add the client with read permissions if they’re copying the data from us or read/write if they’re sending us data - click Send notification to send them the link (If client hasn’t created their Globus signin yet, follow suggestions below for getting the login, then return here to add them)
(Note, sometimes globus doesn’t show the person (ugh!) - try name (first, last, or first last), email address, scrolling, etc and make sure you’re spelling everything correctly. )
6b.Use the Roles tab to add Shannan with administrator permissions so there’s backup if needed
7. Done! Go back to the permissions tab, click on “show link for sharing” copy to clipboard. Add the link to instructions/invites/note it a basecamp message - all helps with followup/followthrough/reminder
8. Suggestions below if the client needs instructions for using Globus.

For more information on creating a globus share:

https://docs.globus.org/guides/tutorials/manage-files/share-files/ 

[Globus Sharing for Storage Providers Youtube](https://www.youtube.com/watch?v=_AkfD881P5A&t=982s)


## Client Globus Connect instructions if needed.

### The first step is telling the client about globus, if needed, and ask them for their globus login or to set it up. Here's a script:

Globus LogIN help 

Go to [https://www.globus.org/](https://www.globus.org/) and click on "LogIn/Use your existing organizational login/Harvard University", this will take you to Harvard Key. If you don’t have Harvard Key you can create a globus id - [https://www.globusid.org/create](https://www.globus.org/create) 
If the client is new to Globus, you can edit the following to give them some instructions. (more examples below)

### Second - once you get their Globus login/email

#### Option 1
First, (ask or check yourself) - some organizations support Globus as HMS does. (ie. BU, Brown, NIH) This makes the transfer easy. Create the guest collection & send the client the link. When they click on the link you provide, our data shows on one side in the globus file manager screen and they enter their organizational info on the other side of the screen. Press start under our guest collection and the data will copy to them. If their org doesn’t have Globus server status, they’ll have to install/use Globus Connect. (Option 2) (https://www.globus.org/globus-connect)

#### Option 2
Following up on the globus notification, here are some instructions for copying the data to your lab. Am I right that you’re copying data to a drive on your computer? (external hard drive, dropbox, etc.-anything mounted is fine. However, we always recommend a server that’s backed up for best data management practice.)
 
1. Mount the drive on your computer if not there already.
2. Install the globus connect app on your system. (https://www.globus.org/globus-connect-personal).
3. Once installed, add/identify the drive where you’re copying the data.

   -3a.Instructions on a mac:  – choose Preferences/Options - use + sign to add destination - this brings up the finder for choosing the destination of your data.

   -3b. On windows: Right-click the Globus Connect Personal icon in the taskbar and select "Options…​". Under the "Access" tab, click the "+" icon and select the folder you wish to make accessible. 
5. In the File Manager window (give link) you’ll see our/HBC “source” (“20250715_PIname_hbc#_”) on the left side. On the other side, choose your destination drive/location using the “search” – then click on whatever you setup using the + above in 3.a.  Click on the START button on the side of our data to copy it to your destination.

I’m happy to get on a zoom call if that helps.


## More sample Globus scripts:
**Brown University example**

For data transfer we use Globus.  

HMS Globus secure transfer service: It works well to transfer the data from O2 to external servers and individual laptops/computers. We need you to sign in and give us your Globus user ID or the email you used to sign in. We'd then share the folder with you to access your data and you will receive an email with instructions.  

It seems Brown uses Globus too, so for logging in: Go to https://www.globus.org/ and click on "LogIn/Use your existing organizational login ". 

You can also download the globus personal app to your computer/laptop but I’d recommend having all the data on a server (with backup). 

https://www.globus.org/globus-connect-personal 

(Additional instructions: How-to: https://docs.globus.org/how-to/get-started/) 

Please, let us know if you have any comments or questions. We can also work with you and/or the IT team if that helps. 

Best,

**Mass General, MGB example**

- not a globus server site so need to use globus connect:
  
1. Mount your ERISTwo drive on your computer.

2. Install the globus connect app on your system. (https://www.globus.org/globus-connect-personal).

3. Once installed, add/identify the ERISTwo endpoint where you’re transferring the data from.

    3a. Instructions on a mac:  – choose Preferences/Options - use + sign to add destination - this brings up the finder for choosing the source of your data.

4. In the File Manager window (give the link)- you’ll see our/HBC “destination” (“PIname_data_HBCcode”) on the left side. On the other side, choose your data source location using the “search” – then click on whatever you setup using the + above in 3.a.  Click on the START button on the side of your data to transfer it to us.

## Misc Globus tips, notes and Troubleshooting:

-keep it simple - ask for login first. When you get that, then start with globus details if needed. 

-In file manager, if needed, the top right, choose the two panel option (middle) to show the source and destination of the transfer. (if your link ends in “path=%2F” it’s all set!-should be the default)) It doesn’t matter which side is the source or destination - just take care to press the Start> button under the right “source”

-cancel the transfer and restart it with the sync option. In the “Transfer and timer options” menu between the two start buttons. You can choose the sync option when you first begin the transfer in the file manager. This will restart the transfer and skip the files that are already transferred. 

-Sorry this took 7 months and ten minutes to complete! – lol. i don’t have an IT specialist in my lab right now and so was a little fearful.  Now i am stronger. 

-for best practices we prefer to hear you’re copying it to a (CHB, BWH, BIDMC, etc.) high performance computing site that’s backed up regularly, etc.  

-Globus Connect: By default, the only folder listed is your home directory. (Ie C:\users\documents) so after clicking the + sign, you may need to go up a level or two to have access to other volumes available. (always leave the home directory.)
