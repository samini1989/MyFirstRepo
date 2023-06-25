namespace DNmodel
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.Runbtn = new System.Windows.Forms.Button();
            this.RunGAbtn = new System.Windows.Forms.Button();
            this.SuspendLayout();
            // 
            // Runbtn
            // 
            this.Runbtn.Location = new System.Drawing.Point(97, 162);
            this.Runbtn.Name = "Runbtn";
            this.Runbtn.Size = new System.Drawing.Size(75, 23);
            this.Runbtn.TabIndex = 1;
            this.Runbtn.Text = "Run";
            this.Runbtn.UseVisualStyleBackColor = true;
            this.Runbtn.Click += new System.EventHandler(this.Runbtn_Click);
            // 
            // RunGAbtn
            // 
            this.RunGAbtn.Location = new System.Drawing.Point(240, 162);
            this.RunGAbtn.Name = "RunGAbtn";
            this.RunGAbtn.Size = new System.Drawing.Size(75, 23);
            this.RunGAbtn.TabIndex = 2;
            this.RunGAbtn.Text = "Run GA";
            this.RunGAbtn.UseVisualStyleBackColor = true;
            this.RunGAbtn.Click += new System.EventHandler(this.RunGAbtn_Click);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(446, 283);
            this.Controls.Add(this.RunGAbtn);
            this.Controls.Add(this.Runbtn);
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);

        }

        #endregion
        private System.Windows.Forms.Button Runbtn;
        private System.Windows.Forms.Button RunGAbtn;
    }
}

